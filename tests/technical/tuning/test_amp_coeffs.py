# Coherent amplitude coefficients through the actual GRANIITTI steering and production cards
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyjson5
import pytest
import torch
from core import resource
from core.io.serialize import load_json_file
from core.numerics.interp import ChebyshevGrid
from core.tune.drivers.graniitti.ampfit.coefficients import (
    ContinuumCoefficients,
    ResonanceCoefficients,
    transfer_weights,
)
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.drivers.graniitti.tunesetup import card
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.parameters.space import normalize_param_space
from core.tune.tunesetup import load_tunesetup

from submit import campaign_source

ROOT = Path(__file__).resolve().parents[3]


# Preserve complex polynomial residuals, beam exchange and derivatives after transfer extraction
@pytest.mark.parametrize("inverse", [False, True])
def test_transfer_interp_complex_residual(inverse):
    axis = ChebyshevGrid(0.2, 5.0, 5)
    transfer = torch.tensor([[-0.04, -0.12], [-0.49, -0.64]], dtype=torch.float64)
    power = 1.3
    poles = axis.nodes if inverse else axis.nodes.reciprocal()
    residual = 1 + 0.2 * axis.nodes + 0.03j * axis.nodes.square()
    amplitudes = (1 - poles[None, :, None] * transfer[:, None, :]).pow(-power).prod(-1) * residual

    # Contract the actual interpolation weights with complex amplitudes on the grid
    def evaluate(value):
        weights = transfer_weights(axis, value, transfer, power, inverse)
        return (weights * amplitudes).sum(-1)

    for point in [*axis.nodes.tolist(), 1.7, 3.8]:
        value = axis.nodes.new_tensor(point)
        pole = value if inverse else value.reciprocal()
        expected = (1 - pole * transfer).pow(-power).prod(-1) * (1 + 0.2 * value + 0.03j * value.square())
        torch.testing.assert_close(evaluate(value), expected, atol=1e-12, rtol=1e-12)
        torch.testing.assert_close(transfer_weights(axis, value, transfer, power, inverse),
                                   transfer_weights(axis, value, transfer.flip(-1), power, inverse))
    value = axis.nodes.new_tensor(1.7, requires_grad=True)
    assert torch.autograd.gradcheck(lambda point: torch.view_as_real(evaluate(point)), (value,))
    assert torch.autograd.gradgradcheck(lambda point: torch.view_as_real(evaluate(point)), (value,))


# Construct the physical pure projector independently from compact magnetic projection rows
def spin_projector(rows, spin):
    state = np.zeros(2 * spin + 1, dtype=complex)
    for m, magnitude, angle in rows:
        state[spin + m] = magnitude * np.exp(1j * angle)
        state[spin - m] = (-1)**m * state[spin + m]
    return np.outer(state, state.conj()) / np.vdot(state, state).real


# Reject mixed spin ensembles before preparing a coherent amplitude bank
def test_ampfit_rejects_density_production(tmp_path):
    driver = GraniittiDriver()
    tune = tmp_path / "density"
    driver.create_steering_card(
        param_space={"RES|f2_1270:MP:polarization.mode": "rho"}, tunename=str(tune), cdir=str(ROOT),
    )
    with pytest.raises(ValueError, match="rho requires a sampling optimizer"):
        ResonanceCoefficients(driver=driver, path=tune, model="MP", resonances=["f2_1270"],
                              pid=[211, -211], initial={})


# Reject pole variations whose decay or line shape cannot factor out of the native amplitudes
@pytest.mark.parametrize("parameters,reason", [
    ({"RES|f2_1950:GP:BW": "running-width"}, "supported line shape"),
    ({"REGGE|DECAY_BARRIERS.GP": True}, "decay barriers"),
])
def test_pole_requires_supported_decay(tmp_path, parameters, reason):
    driver = GraniittiDriver()
    driver.create_steering_card(param_space=parameters, tunename=str(tmp_path), cdir=str(ROOT))
    resonance = load_json_file(tmp_path / "RES/f2_1950.json", loader=pyjson5.load)["PARAM_RES"]
    with pytest.raises(ValueError, match=reason):
        ResonanceCoefficients(driver=driver, path=tmp_path, model="GP", resonances=["f2_1950"], pid=[211, -211],
                              initial={"RES|f2_1950:GP:width": resonance["MODELS"]["GP"]["width"]},
                              decay=load_settings(resource("tune/settings/ampfit.json"))["bank"])


# Prepare the source tune and reject unsupported active coordinates before generator jobs
@pytest.mark.parametrize("model", ["MP", "XP", "GP"])
def test_preflight_plans_coefficients(tmp_path, model):
    driver = GraniittiDriver()
    dataset = "icepack/SOFTCEP/CMS_2752118/pipi_0p7/dataset.json"
    key = f"RES|f0_500:{model}:phi"
    space = {key: {"type": "uniform", "lower": -math.pi, "upper": math.pi}}
    tunesetup = SimpleNamespace(
        param_space=space, aux_param_space={},
        datacards=[card(dataset, nevents=1, loopscreen=True, xsmode="sample", swap_process=f"{model}[RES+CON]<F>")],
    )
    args = SimpleNamespace(algorithm="ampfit", cost="chi2", cdir=str(ROOT), tune_default="TUNE0",
                           run_name=str(tmp_path / "preflight"), tunesetup="amplitude_preflight", obs_module="default",
                           data_covariance_mode="diagonal", mc_correlation_events=None, mc_correlation_weighting="sample",
                           ampfit_settings=load_settings(resource("tune/settings/ampfit.json")))
    driver.prepare_tunesetup(tunesetup=tunesetup, args=args)
    source = tmp_path / "preflight" / "init" / "ampfit_preflight"
    assert source.is_dir()
    expected = driver.get_initial_param(space, {}, str(ROOT))[key]
    actual = load_json_file(source / "RES/f0_500.json")["PARAM_RES"]["MODELS"][model]["phi"]
    assert actual == pytest.approx(expected)

    tunesetup.param_space["REGGE|s0"] = {"type": "uniform", "lower": 0.1, "upper": 10.0}
    with pytest.raises(ValueError, match="not represented in the amplitude bank.*s0"):
        driver.prepare_tunesetup(tunesetup=tunesetup, args=args)


# Compare complex production couplings with the normal card writer for every fitted production model
@pytest.mark.parametrize("model,tunesetup", [("MP", "tune-mpom-res-con-coh"), ("XP", "tune-xpom-res-con"), ("GP", "tune-gpom-res-con")])
def test_res_coeffs_match_steering(tmp_path, model, tunesetup):
    driver = GraniittiDriver()
    tune = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(tunesetup))
    initial = driver.get_initial_param(tune.param_space, tune.aux_param_space, str(ROOT))
    names = sorted({key.split("|")[1].split(":")[0] for key in initial if key.startswith("RES|")})
    coefficients = ResonanceCoefficients(driver=driver, path=ROOT / "modeldata/TUNE0", model=model,
                                        resonances=names, pid=[211, -211], initial=initial | tune.aux_param_space,
                                        decay=load_settings(resource("tune/settings/ampfit.json"))["bank"])
    target = tmp_path / model
    driver.create_steering_card(param_space=initial | tune.aux_param_space, tunename=str(target), cdir=str(ROOT))
    mass2 = torch.tensor([coefficients.resonances[name]["card"]["PARAM_RES"]["MODELS"][coefficients.model]["mass"]**2 for name in names], dtype=torch.float64)
    # At the source pole the line shape ratio is unity, independent of the daughter masses
    values = coefficients.evaluate(initial | tune.aux_param_space, mass2, [0.0, 0.0], torch.zeros((len(mass2), 2), dtype=mass2.dtype)).detach().numpy()
    for index, column in enumerate(coefficients.columns):
        name = column["resonance"]
        card = load_json_file(target / "RES" / f"{name}.json", loader=pyjson5.load)
        res = card["PARAM_RES"]
        block = driver._resonance_block(card, model)
        expected = np.exp(1j * res["MODELS"][model]["phi"])
        if column["magnitude"] is not None:
            field = column["magnitude"].rsplit(":", 1)[1]
            label, indices = driver._indexed_parameter(field)
            coupling = block[label]
            for row in indices[:-1]:
                coupling = coupling[row]
            position = indices[-1]
            expected *= coupling[position] * np.exp(1j * coupling[position + 1])
        if "spin" not in column:
            np.testing.assert_allclose(values[names.index(name), index], expected, atol=1e-11)
            continue
        spin = res["spinX2"] // 2
        actual = np.zeros((2 * spin + 1, 2 * spin + 1), dtype=complex)
        for position, component in enumerate(coefficients.columns):
            if component["resonance"] == name:
                rows = component["parameters"][f"RES|{name}:MP:polarization.a_Jz"]
                actual += values[names.index(name), position] * spin_projector(rows, spin)
        np.testing.assert_allclose(actual, expected * spin_projector(block["polarization"]["a_Jz"], spin), atol=1e-11)


# Differentiate coupled production, pole and form factor coordinates through one complete coefficient map
@pytest.mark.parametrize("model,tunesetup", [("MP", "tune-mpom-res-con-coh"), ("XP", "tune-xpom-res-con"), ("GP", "tune-gpom-res-con")])
def test_resonance_coefficients_autograd(model, tunesetup):
    driver = GraniittiDriver()
    tune = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(tunesetup))
    initial = driver.get_initial_param(tune.param_space, tune.aux_param_space, str(ROOT)) | tune.aux_param_space
    initial = {key: value for key, value in initial.items() if key.startswith("RES|f0_500:")}
    coefficients = ResonanceCoefficients(driver=driver, path=ROOT / "modeldata/TUNE0", model=model,
                                        resonances=["f0_500"], pid=[211, -211], initial=initial,
                                        decay=load_settings(resource("tune/settings/ampfit.json"))["bank"])
    keys = sorted(key for key, value in initial.items() if isinstance(value, float))
    theta = torch.tensor([initial[key] for key in keys], dtype=torch.float64, requires_grad=True)
    mass2 = theta.new_tensor([0.3, 0.7, 1.2])

    # Test the full complex coefficient map without substituting an analytic Jacobian
    def evaluate(value):
        parameters = initial | dict(zip(keys, value, strict=True))
        return torch.view_as_real(coefficients.evaluate(parameters, mass2, [0.0, 0.0], torch.zeros((len(mass2), 2), dtype=mass2.dtype)))

    assert torch.autograd.gradcheck(evaluate, (theta,), eps=1e-5)
    assert torch.autograd.gradgradcheck(evaluate, (theta,), eps=1e-5)


# Differentiate complex continuum coefficients at interior points of the fitted ranges
@pytest.mark.parametrize("model,tunesetup", [("MP", "tune-mpom-res-con-coh"), ("XP", "tune-xpom-res-con"), ("GP", "tune-gpom-res-con")])
@pytest.mark.parametrize("pid", [211, 321])
def test_continuum_coefficients_autograd(model, tunesetup, pid):
    driver = GraniittiDriver()
    tune = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(tunesetup))
    original = driver.get_initial_param(tune.param_space, tune.aux_param_space, str(ROOT))
    bounds = normalize_param_space(tune.param_space)
    initial = {key: value for key, value in (original | tune.aux_param_space).items()
               if key.startswith(f"CON_{model}|") and f":[{pid},{pid}]" in key}
    coefficients = ContinuumCoefficients(driver=driver, path=ROOT / "modeldata/TUNE0", model=model,
                                         pid=[pid, -pid], initial=initial, bounds=bounds, nodes=4)
    keys = sorted(key for key in initial if key in bounds)
    theta = torch.tensor([(bounds[key]["lower"] + bounds[key]["upper"]) / 2 for key in keys],
                         dtype=torch.float64, requires_grad=True)
    mass2 = theta.new_tensor([0.8, 2.0, 5.0])

    # Differentiate exchange products, spline weights and the common central veto together
    def evaluate(value):
        return torch.view_as_real(coefficients.evaluate(initial | dict(zip(keys, value, strict=True)), mass2, mass2, -mass2[:, None].expand(-1, 2)).to(torch.complex128))

    assert torch.autograd.gradcheck(evaluate, (theta,), eps=1e-5, fast_mode=True)
    assert torch.autograd.gradgradcheck(evaluate, (theta,), eps=1e-5, fast_mode=True)


# Differentiate the supported meson form factors through actual source card steering
@pytest.mark.parametrize("ff_type", ["exp", "power", "orear", "logexp"])
@pytest.mark.parametrize("pid", [211, 321])
def test_continuum_form_factor_autograd(tmp_path, ff_type, pid):
    driver = GraniittiDriver()
    prefix = f"CON_GP|990:[{pid},{pid}]:"
    form = driver._continuum_ff_template(ff_type)
    keys = [prefix + "FF_offshell." + field for field in form if field not in {"type", "norm"}]
    keys.append(prefix + "FF_transfer.LambdaInv2")
    fixed = {prefix + "FF_offshell.type": ff_type}
    initial = driver.get_initial_param(dict.fromkeys(keys), fixed, str(ROOT)) | fixed
    driver.create_steering_card(param_space=initial, tunename=str(tmp_path), cdir=str(ROOT))
    bounds = {key: {"lower": initial[key] * 0.9, "upper": initial[key] * 1.1} for key in keys}
    coefficients = ContinuumCoefficients(driver=driver, path=tmp_path, model="GP", pid=[pid, -pid],
                                         initial=initial, bounds=bounds, nodes=3)
    theta = torch.tensor([initial[key] * 1.04 for key in keys], dtype=torch.float64, requires_grad=True)
    mass2 = theta.new_tensor([2.0, 5.0])

    # Differentiate the complex weights away from interpolation nodes
    def evaluate(value):
        weights = coefficients.evaluate(initial | dict(zip(keys, value, strict=True)), mass2, mass2, -mass2[:, None].expand(-1, 2))
        return torch.view_as_real(weights.to(torch.complex128))

    assert torch.autograd.gradcheck(evaluate, (theta,), fast_mode=True)
    assert torch.autograd.gradgradcheck(evaluate, (theta,), fast_mode=True)


# Differentiate direct spin coordinates while the native amplitude supplies the diphoton width
# The unit spin components must retain the same physical production normalization
def test_diphoton_polarization_coefficients(tmp_path):
    from core.tune.parameters.tools import spherical_angle_key

    driver = GraniittiDriver()
    name = "f2_1270_yy"
    prefix = f"RES|{name}:MP:"
    angle = spherical_angle_key(prefix + "polarization.a_Jz", 0)
    initial = {angle: 0.43, prefix + "g[1]": 0.21, "REGGE|MP_FRAME": "CS"}
    coefficients = ResonanceCoefficients(driver=driver, path=ROOT / "modeldata/TUNE0", model="MP",
                                        resonances=[name], pid=[211, -211], initial=initial)
    mass2 = torch.tensor([coefficients.resonances[name]["card"]["PARAM_RES"]["MODELS"][coefficients.model]["mass"]**2], dtype=torch.float64)
    theta = torch.tensor([initial[angle], initial[prefix + "g[1]"]], dtype=torch.float64, requires_grad=True)

    # Evaluate the production state using the same coefficient map as ampfit
    def evaluate(value):
        return coefficients.evaluate(initial | {angle: value[0], prefix + "g[1]": value[1]}, mass2, [0.0, 0.0], torch.zeros((len(mass2), 2), dtype=mass2.dtype))

    weights = evaluate(theta).detach().numpy()[0] * np.exp(-1j * initial[prefix + "g[1]"])
    density = sum(weight * spin_projector(column["parameters"][prefix + "polarization.a_Jz"], 2)
                  for weight, column in zip(weights, coefficients.columns, strict=True))
    np.testing.assert_allclose(density @ density, density, atol=1e-12)
    assert np.trace(density).real == pytest.approx(1.0)
    assert torch.autograd.gradcheck(lambda value: torch.view_as_real(evaluate(value)), (theta,))
    for index, column in enumerate(coefficients.columns):
        target = tmp_path / f"diphoton_{index}"
        driver.create_steering_card(param_space={"REGGE|MP_FRAME": "CS"} | column["parameters"], tunename=str(target), cdir=str(ROOT))
        prepared = load_json_file(target / "RES" / f"{name}.json")
        assert driver._resonance_block(prepared, "MP")["g"][0] is None


# Compare complex pole reweighting with its analytic ratio for independent model parameters
@pytest.mark.parametrize("model", ["MP", "XP", "GP"])
@pytest.mark.parametrize("mode", ["fixed-width", "kinematic-width"])
def test_model_pole_coefficients(tmp_path, model, mode):
    driver = GraniittiDriver()
    source = load_json_file(ROOT / "modeldata/TUNE0/RES/f0_500.json", loader=pyjson5.load)["PARAM_RES"]
    parameters = {f"RES|f0_500:{name}:{field}": pole[field] * (1.0 + 0.1 * index)
                  for index, (name, pole) in enumerate(source["MODELS"].items()) for field in ("mass", "width")}
    prefix = f"RES|f0_500:{model}:"
    parameters.update({prefix + "BW": mode, f"REGGE|DECAY_BARRIERS.{model}": False,
                       f"DECAY|{source['PDG']}:[211,-211]:FF_decay.{model}": {"type": "none"}})
    target = tmp_path / "poles"
    driver.create_steering_card(param_space=parameters, tunename=str(target), cdir=str(ROOT))
    initial = {key: value for key, value in parameters.items() if key.endswith((":mass", ":width"))}
    coefficients = ResonanceCoefficients(driver=driver, path=target, model=model, resonances=["f0_500"],
                                        pid=[211, -211], initial=initial,
                                        decay=load_settings(resource("tune/settings/ampfit.json"))["bank"])
    assert coefficients.used == {prefix + "mass", prefix + "width"}
    mass0, width0 = initial[prefix + "mass"], initial[prefix + "width"]
    mass, width = 1.03 * mass0, 1.07 * width0
    mass2 = torch.tensor([0.9, 1.1, 1.4], dtype=torch.float64) * mass0**2
    transfer = torch.zeros((len(mass2), 2), dtype=mass2.dtype)
    original = coefficients.evaluate(initial, mass2, [0.0, 0.0], transfer)
    varied = {key: 1.2 * value for key, value in initial.items()}
    varied.update({prefix + "mass": mass0, prefix + "width": width0})
    torch.testing.assert_close(coefficients.evaluate(varied, mass2, [0.0, 0.0], transfer), original)
    varied.update({prefix + "mass": mass, prefix + "width": width})
    imaginary0 = torch.sqrt(mass2) if mode == "kinematic-width" else mass0
    imaginary = torch.sqrt(mass2) if mode == "kinematic-width" else mass
    ratio = (mass2 - mass0**2 + 1j * imaginary0 * width0) / (mass2 - mass**2 + 1j * imaginary * width)
    ratio *= math.sqrt(width * mass / (width0 * mass0))
    form = coefficients.resonances["f0_500"]["form"]
    ratio *= torch.exp(((mass2 - mass0**2)**2 - (mass2 - mass**2)**2) / form["Lambda2"]**2)
    torch.testing.assert_close(coefficients.evaluate(varied, mass2, [0.0, 0.0], transfer), original * ratio[:, None])


# Preserve complex contractions and derivatives when sharing couplings across interpolation nodes
@pytest.mark.parametrize("model", ["MP", "XP", "GP"])
def test_resonance_contraction(model):
    driver = GraniittiDriver()
    resonances = ["f0_500", "f2_1270"]
    keys = [f"RES|{name}:{model}:{field}" for name in resonances for field in ("phi", "FF_transfer.LambdaInv2")]
    initial = driver.get_initial_param(dict.fromkeys(keys), {}, str(ROOT))
    bounds = {key: {"lower": initial[key] * 0.8, "upper": initial[key] * 1.2}
              for key in keys if key.endswith("LambdaInv2")}
    coefficients = ResonanceCoefficients(driver=driver, path=ROOT / "modeldata/TUNE0", model=model,
        resonances=resonances, pid=[211, -211], initial=initial, bounds=bounds, nodes=3)
    theta = torch.tensor(list(initial.values()), dtype=torch.float64, requires_grad=True)
    mass2 = theta.new_tensor([1.0, 2.0, 3.0])
    transfer = -mass2[:, None].expand(-1, 2)
    basis = torch.randn(len(coefficients.columns), len(mass2), 2, dtype=torch.complex128,
                        generator=torch.Generator().manual_seed(917))

    # Contract the actual coefficient map with a general complex external helicity basis
    def evaluate(values, contracted):
        return coefficients.evaluate(dict(zip(initial, values, strict=True)), mass2, [0.0, 0.0], transfer,
                                     basis if contracted else None)

    actual = evaluate(theta, True)
    expected = torch.einsum("cnh,nc->nh", basis, evaluate(theta, False))
    torch.testing.assert_close(actual, expected)
    for part in ("real", "imag"):
        a, = torch.autograd.grad(getattr(actual, part).sum(), theta, retain_graph=True)
        b, = torch.autograd.grad(getattr(expected, part).sum(), theta, retain_graph=True)
        torch.testing.assert_close(a, b)
    assert torch.autograd.gradgradcheck(lambda values: torch.view_as_real(evaluate(values, True)), (theta,), fast_mode=True)
