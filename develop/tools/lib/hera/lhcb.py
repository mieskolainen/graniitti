# Fit HERA and LHCb measurements with direct photoproduction cross sections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from dataclasses import replace
from pathlib import Path

import numpy as np
import pyjson5
from core.io.hepdata_reader import load_dataset_reader
from core.io.serialize import load_json_file
from core.io.steering import load_dataset, resolve_gencard_reference
from core.stats.uncertainty import source_covariance
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import least_squares
from scipy.stats import chi2

from icepack.PHOTOPROD.LHCb_1373746.lhcb_upsilon import energy_fractions

from . import fit as hera_fit
from . import model, pp


# Read the original LHCb spectrum and prepare its direct photon-flux integration
class MeasurementPP:
    """Rapidity spectrum at one or two beam energies"""

    # Resolve all measurement and generator references with the shared icepack functions
    def __init__(self, spec, proton, general, numerics, physics):
        dataset, path = load_dataset(spec['dataset'], cdir=model.REPOSITORY_ROOT)
        panel = dataset['sets'][spec['set']]
        hist = panel['hist'][0]
        if panel['region'] != 'central_system' or not hist.get('differential', True):
            raise ValueError('The direct pp fit requires a differential central-system rapidity measurement')
        self.data = load_dataset_reader(spec['dataset'], cdir=str(model.REPOSITORY_ROOT)).read(
            str(model.REPOSITORY_ROOT / dataset['datapath'] / hist['file']), dataset=panel)
        self.covariance = sum(source_covariance(source) for source in self.data['uncertainties'])
        self.chol = cholesky(self.covariance, lower=True)
        x, self.weights = pp.gauss(numerics['rapidity'], -1, 1)
        low, high = self.data['binedges'].T
        rapidity = ((low[:, None] + high[:, None] + (high - low)[:, None] * x) / 2).ravel()
        self.weights /= 2
        # Predictions describe stable vector mesons, so no decay branching fraction enters
        self.scale = 1e-6 / hist['scale']
        self.spectra, energies = [], []
        tune = Path(proton.source.split(':', 1)[0]).parent
        resonances = [load_json_file(p, loader=pyjson5.load)['PARAM_RES'] for p in (tune / 'RES').glob('*.json')]
        mass = next(r['MODELS']['GP']['mass'] for r in resonances if r['PDG'] == spec['pdg'])
        self.mass, self.rapidity = mass, rapidity
        self.dataset = str(path)
        for name in spec['samples']:
            sample = next(row for row in dataset['samples'] if row['name'] == name)
            card = load_json_file(resolve_gencard_reference(sample['gencard'], dataset_path=path, cdir=model.REPOSITORY_ROOT), loader=pyjson5.load)
            energy = sum(sample['parameters'].get('SCATTERING.ENERGY', card['SCATTERING']['ENERGY']))
            energies.append(energy)
            self.spectra.append(pp.Spectrum(energy, rapidity, mass, proton, general, numerics, physics))
        self.energies = energies
        self.quadrature = numerics
        self.fractions = [1.0] if len(energies) == 1 else energy_fractions(energies)

    # Average differential cross sections in microbarn over the measured beam luminosities
    def differential(self, point, channel, screened=True):
        forward, slope, delta, shrinkage = point[:4]
        curve = (self.mass**2, point[4]) if len(point) > 4 else channel.dlog
        return sum(fraction * spectrum.predict(forward, slope, delta, shrinkage, channel.w0_gev, curve, screened)
                     for fraction, spectrum in zip(self.fractions, self.spectra, strict=True))

    # Average the differential prediction over the measured rapidity bins
    def predict(self, point, channel, screened=True):
        return self.differential(point, channel, screened).reshape(-1, len(self.weights)) @ self.weights * self.scale


# Fit the same absolute gamma-p amplitude to independent HERA and LHCb measurements
# [REFERENCE: arXiv:2409.03496, arXiv:1505.08139]
def fit_channel(channel, proton, measurement, controls):
    forward = channel.forward_dsigma_dt_ub_per_gev2.central
    slope = channel.vector_cross_section_slope_per_gev2
    delta = 4 * (channel.alpha0.central - 1)
    shrinkage = channel.alpha_prime_per_gev2.central
    initial = [forward, slope, delta, shrinkage]
    if channel.pdg == 443:
        hera = hera_fit.jpsi_spectra()
        observed, chol, relative, _ = hera_fit.jpsi_errors(
            [row[1]['y'][row[1]['fit_mask']] for row in hera],
            [row[1]['stat_cov'][np.ix_(row[1]['fit_mask'], row[1]['fit_mask'])] for row in hera],
            [{key: value[row[1]['fit_mask']] for key, value in row[-1].items()} for row in hera])
        initial.append(channel.dlog[1] if channel.dlog else 0.0)
        active = np.arange(5)
    else:
        hera = channel.fit
        observed = np.asarray(hera['data'])
        energy = np.asarray(hera['W'])
        transfer = np.asarray([model.profile_moments(proton, slope + 4 * shrinkage * np.log(w / channel.w0_gev), None)[0] for w in energy])
        if channel.dlog is not None:
            transfer *= model.dlog_factor(energy, channel.w0_gev, *channel.dlog)
        active = np.array([0, 2])
    initial = np.asarray(initial)

    # Keep HERA nuisance constraints and the full measured LHCb covariance
    def residual(point):
        pars = initial.copy()
        pars[active] = point
        if channel.pdg == 443:
            gamma_p = np.concatenate([hera_fit.jpsi_spectrum([pars[0] * 1000, *pars[1:4]], proton, d, w, flux,
                                     (measurement.mass**2, pars[4]))[d['fit_mask']] for _, d, w, flux, *_ in hera])
            gamma_p = hera_fit.profiled_residual(gamma_p, observed, chol, relative)
        else:
            difference = pars[0] * (energy / channel.w0_gev)**pars[2] * transfer - observed
            gamma_p = difference / np.where(difference >= 0, hera['errors_up'], hera['errors_down'])
        pp = measurement.predict(pars, channel)
        return np.r_[gamma_p, solve_triangular(measurement.chol, pp - measurement.data['y'], lower=True)]

    upper = np.full(len(active), np.inf)
    if channel.pdg == 443:
        upper[-1] = controls['curvature_max']
    result = least_squares(residual, initial[active], bounds=(0, upper), x_scale='jac',
                           **{key: controls[key] for key in ('ftol', 'xtol', 'gtol', 'max_nfev')})
    if not result.success or np.linalg.matrix_rank(result.jac) != len(active):
        raise ValueError('Direct HERA/LHCb fit failed or is not identifiable: ' + result.message)
    point = initial.copy()
    point[active] = result.x
    cov = np.zeros((len(point), len(point)))
    cov[np.ix_(active, active)] = np.linalg.inv(result.jac.T @ result.jac)
    error = np.sqrt(np.diag(cov))
    values = residual(result.x)
    prediction = measurement.predict(point, channel)
    assumptions = [*channel.assumptions,
                   'Direct bin integration of the absolute gamma-p amplitude and transverse Dirac-Pauli photon current',
                   'Small-x on-shell vectors, neglecting photon virtuality and nonlogarithmic EPA terms',
                   'Single-channel SOFT absorption with coherent interference between the two photon directions',
                   'No generated events, amplitude banks, free pp normalization or fitted survival factor',
                   'Covariance conditional on the fixed SOFT model and the stated high-energy approximation']
    if channel.pdg == 553:
        assumptions.append('Upsilon transfer profile fixed by HERA, common beam efficiencies in the combined sample')
    ndf = len(observed) + len(prediction) - len(active)
    report = {'pvalue': float(chi2.sf(values @ values, ndf)), 'parameters': point.tolist(), 'parameter_order': ['forward_ub_per_GeV2', 'B_vector', 'delta_forward', 'alpha_prime']
              + (['double_log_c'] if len(point) == 5 else []), 'covariance': cov.tolist(), 'chi2': float(values @ values),
              'ndf': ndf, 'hera_chi2': float(values[:-len(prediction)] @ values[:-len(prediction)]),
              'lhcb_chi2': float(values[-len(prediction):] @ values[-len(prediction):]), 'assumptions': assumptions,
              'quadrature': measurement.quadrature,
              'lhcb': {'dataset': measurement.dataset, 'energies_GeV': measurement.energies,
                       'mass_GeV': measurement.mass, 'rapidity': measurement.rapidity.tolist(),
                       'curve': (measurement.differential(point, channel) * measurement.scale).tolist(),
                       'fractions': list(measurement.fractions), 'binedges': measurement.data['binedges'].tolist(), 'data': measurement.data['y'].tolist(),
                       'prediction': prediction.tolist(), 'born': measurement.predict(point, channel, False).tolist(),
                       'errors': np.sqrt(np.diag(measurement.covariance)).tolist()}, 'evaluations': result.nfev, 'active_bounds': result.active_mask.tolist()}
    if channel.pdg == 443:
        report['spectra'] = [{'spectrum': name, 'source': d['source'], 'binedges': d['binedges'].tolist(), 'data': d['y'].tolist(),
                             'prediction': hera_fit.jpsi_spectrum([point[0] * 1000, *point[1:4]], proton, d, w, flux,
                                                         (measurement.mass**2, point[4])).tolist(),
                             'errors': np.sqrt(np.diag(c)).tolist(), 'shape_fit_bins': d['fit_mask'].tolist()}
                            for name, d, w, flux, c, _ in hera]
    else:
        report['energy'] = {'W': energy.tolist(), 'data': observed.tolist(),
                            'prediction': (point[0] * (energy / channel.w0_gev)**point[2] * transfer).tolist(),
                            'errors': ((np.asarray(hera['errors_up']) + hera['errors_down']) / 2).tolist()}
    return replace(channel, forward_dsigma_dt_ub_per_gev2=model.Measurement(point[0], error[0], error[0], 0, 0, 'fit'),
                   vector_cross_section_slope_per_gev2=point[1], alpha0=model.Measurement(1 + point[2] / 4, error[2] / 4, error[2] / 4, 0, 0, 'fit'),
                   alpha_prime_per_gev2=model.Measurement(point[3], error[3], error[3], 0, 0, 'fit' if channel.pdg == 443 else 'fixed'),
                   normalization_origin='Joint absolute HERA and LHCb direct fit', trajectory_origin='HERA and LHCb W dependence',
                   fit=report, assumptions=tuple(assumptions), dlog=(measurement.mass**2, point[4]) if len(point) == 5 else channel.dlog)


# Derive the heavy-vector channels and retain explicit excited-state assumptions
def fit(channels, proton, config):
    general = load_json_file(proton.source.split(':', 1)[0], loader=pyjson5.load)
    physics = config['physics']
    for spec in config['channels']:
        parent = channels[spec['pdg']]
        measurement = MeasurementPP(spec, proton, general, config['quadrature'], physics)
        fitted = fit_channel(parent, proton, measurement, config['fit'])
        channels[spec['pdg']] = fitted
        if spec['pdg'] == 553:
            for pdg in (100553, 200553):
                old = channels[pdg]
                ratio = old.forward_dsigma_dt_ub_per_gev2.central / parent.forward_dsigma_dt_ub_per_gev2.central
                channels[pdg] = replace(fitted, key=old.key, name=old.name, pdg=pdg, fit={'reference_pdg': 553},
                                       forward_dsigma_dt_ub_per_gev2=model.rescale(fitted.forward_dsigma_dt_ub_per_gev2, ratio))
