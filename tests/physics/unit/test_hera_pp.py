# Check direct HERA pp quadrature, complex screening phases and beam symmetry
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
from dataclasses import replace
from pathlib import Path

import numpy as np
import pyjson5
import pytest
from core.io.serialize import load_json_file
from scipy.constants import physical_constants

from develop.tools.lib import SETTINGS
from develop.tools.lib.hera.lhcb import MeasurementPP
from develop.tools.lib.hera.model import load_proton, proton_profile
from develop.tools.lib.hera.pp import Spectrum, electromagnetic, masses, smatrix

ROOT = Path(__file__).resolve().parents[3]


# Read model, quadrature and data settings from their actual repository cards
@pytest.fixture(scope='module')
def inputs():
    path = ROOT / 'modeldata/TUNE0/GENERAL.json'
    config = load_json_file(SETTINGS)['hera']
    return load_proton(path), load_json_file(path, loader=pyjson5.load), config, config['physics']


# Check complex Fourier normalization and rotating phases against the Gaussian transform
@pytest.mark.parametrize('unitarization', ['exp', 'q_exp'])
def test_soft_gaussian(inputs, unitarization):
    proton, source, config, _ = inputs
    general = copy.deepcopy(source)
    settings = general['PARAM_SOFT']
    settings['EXCHANGE_DEF']['P']['trajectory_mode'] = 'linear'
    model = settings['MODEL'][settings['active_model']]
    model['EIKONAL']['unitarization'] = unitarization
    slope = 2 * proton.amplitude_slope_per_gev2
    proton = replace(proton, form_factor='EXP', parameters=((slope,),))
    alpha0, ap = model['EXCHANGE']['P']['alpha']
    energy = next(row[1] for row in general['PARAM_REGGE']['photoprod'] if row[0] == 443)
    b = np.linspace(0, config['quadrature']['impact_max'], 31)
    width = slope + 2 * ap * np.log(energy) - 0.5j * np.pi * ap
    coefficient = (-model['EXCHANGE']['P']['sign'] * np.exp(-0.5j * np.pi * alpha0)
                   * energy**(2 * (alpha0 - 1)) * proton.beam_residue_per_gev**2)
    chi = coefficient / (8 * np.pi * width) * np.exp(-b*b / (4 * width))
    q = model['EIKONAL']['q']
    expected = np.exp(1j * chi) if unitarization == 'exp' else (1 + (1-q)*1j*chi)**(1/(1-q))
    np.testing.assert_allclose(smatrix(b, energy, proton, general, config['quadrature']), expected, atol=1e-11)


# Check proton charge and anomalous magnetic moment at zero photon virtuality
@pytest.mark.parametrize('form', ['DIPOLE', 'KELLY'])
def test_photon_form(inputs, form):
    *_, physics = inputs
    f1, f2 = electromagnetic(0, masses()[2212], form, physics)
    assert f1 == pytest.approx(1)
    assert f2 == pytest.approx(physical_constants['proton mag. mom. to nuclear magneton ratio'][0] - 1)


# Exercise the real icepack reader, unit conversion, quadrature and pp photon directions
@pytest.mark.parametrize('index', [0, 1])
def test_pp_spectrum(inputs, index):
    proton, general, config, physics = inputs
    spec = config['channels'][index]
    measurement = MeasurementPP(spec, proton, general, config['quadrature'], physics)
    row = next(row for row in general['PARAM_REGGE']['photoprod'] if row[0] == spec['pdg'])
    _, w0, slope, alpha0, ap = row
    spectrum = measurement.spectra[0]
    y = np.log(spectrum.x[0, len(spectrum.x[0]) // 2] * spectrum.energy / measurement.mass)
    rapidity = np.array([-y, 0, y])
    direct = Spectrum(spectrum.energy, rapidity, measurement.mass, proton, general, config['quadrature'], physics)
    born = direct.predict(1, slope, 4*(alpha0-1), ap, w0, screened=False)
    screened = direct.predict(1, slope, 4*(alpha0-1), ap, w0)
    assert np.all((screened > 0) & (screened < born))
    np.testing.assert_allclose(screened, screened[::-1], rtol=1e-12)
    np.testing.assert_allclose(direct.predict(2, slope, 4*(alpha0-1), ap, w0), 2*screened, rtol=1e-12)
    finer = dict(config['quadrature'])
    for key in ('momentum', 'impact', 'angle', 'soft_nodes'):
        finer[key] *= 2
    finer['impact_max'] *= 1.5
    refined = Spectrum(spectrum.energy, rapidity, measurement.mass, proton, general, finer, physics)
    np.testing.assert_allclose(refined.predict(1, slope, 4*(alpha0-1), ap, w0), screened, rtol=2e-5)
    # Parseval relates the complex hard amplitude in transfer and impact parameter
    profile = proton_profile(direct.q**2, proton) * np.exp(-0.5*(slope - 1j*np.pi*ap)*direct.q**2)
    hard = profile @ direct.hankel
    np.testing.assert_allclose(2*np.sum(direct.rw*abs(hard)**2),
                               2*np.sum(direct.q*direct.qw*abs(profile)**2), rtol=1e-9)
    free = Spectrum(spectrum.energy, rapidity, measurement.mass, replace(proton, beam_residue_per_gev=0),
                    general, config['quadrature'], physics)
    np.testing.assert_allclose(free.predict(1, slope, 4*(alpha0-1), ap, w0), born, rtol=1e-12)
