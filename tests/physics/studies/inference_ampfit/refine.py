# Refine screened GP amplitudes against native STAR and CMS histogram comparisons
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import copy
import json
import pickle
import time
from pathlib import Path

import numpy as np
import torch
from core.stats import hist as histogram
from core.stats import objective
from core.tune.drivers.graniitti.ampfit.amplitude import AmplitudeBank
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.tunesetup import load_tunesetup
from scipy.optimize import minimize

from submit import campaign_source


# Load the prepared screened amplitudes with the normal icepack cuts and observables
def load_banks(root, run, controls):
    tunesetup = load_tunesetup(cdir=root, simdriver='graniitti',
                          name=campaign_source("tune-gpom-star-cms-ampfit"))
    driver = GraniittiDriver()
    driver.init_data(str(run), tunesetup.datacards, 'default', str(root), pickle_dump=False)
    banks = []
    for index in (0, 1, 2):
        indices = driver._sample_set_indices(index=index, sample_name='nominal')
        banks.append(AmplitudeBank(driver=driver, directory=run / 'results/amplitude' / str(index) / 'nominal',
            obs=[driver.obs[index][i] for i in indices], pid=[driver.pid[index][i] for i in indices],
            cuts=[driver.cuts[index][i] for i in indices], controls=controls))
    return banks


# Contract the fixed resonance amplitudes once while retaining continuum interference
class PionResponse:
    # Prepare the exact common resonance contribution at the fitted point
    def __init__(self, bank, parameters):
        self.bank = bank
        with torch.no_grad():
            coefficients = bank.coefficients(parameters)
            size = len(bank.resonance.columns)
            coefficients[:, size:] = 0
            self.resonances = bank.coherent(coefficients)
        self.continuum = bank.amplitudes[1 + size:]

    # Compute native histogram changes from coherent continuum variations
    def predict(self, parameters):
        bank = self.bank
        coefficients = bank.continuum.evaluate(parameters, bank.mass2, bank.virtuality, bank.transfer)
        parts = []
        for start in range(0, len(bank.mass2), bank.batch_events):
            stop = start + bank.batch_events
            parts.append((self.continuum[:, start:stop] * coefficients[start:stop].T[..., None]).sum(0))
        amplitude = self.resonances + torch.cat(parts)
        ratio = amplitude.abs().square().sum(1) / bank.reference
        return bank.weighted_histograms(ratio)


# Add the common-event amplitude response to the independent native prediction
def corrected_histograms(prediction, anchor, native):
    result = []
    for current, source, reference in zip(prediction, anchor, native, strict=True):
        histograms = {}
        for name, item in current.items():
            hist, old, actual = item['hdata'], source[name]['hdata'], reference[name]['hdata']
            counts = hist.counts_scaled - old.counts_scaled + torch.as_tensor(actual.counts_scaled)
            histograms[name] = {**item, 'hdata': histogram.hobj(counts=counts, errs=torch.as_tensor(actual.errs_scaled),
                bins=actual.bins, cbins=actual.cbins, valid=actual.valid)}
        result.append(histograms)
    return result


# Sum ordinary icepack chi2 values and isolate the inclusive STAR pion mass tail
def residuals(mc, data, tail_min):
    total, bins = torch.tensor(0.0, dtype=torch.float64), 0
    for prediction, measured in zip(mc, data, strict=True):
        for name, entry in measured.items():
            if entry['fitw'] <= 0:
                continue
            chi2, ndf = objective.chi2_cost(prediction[name]['hdata'], entry['hdata'])
            total = total + chi2
            bins += ndf
    if bins < 1:
        raise FloatingPointError('No finite histogram residuals at the proposed amplitude point')
    if 'M' not in data[0]:
        return total / bins, total * 0
    hist, measured = mc[0]['M']['hdata'], data[0]['M']['hdata']
    selected = np.asarray(hist.valid & measured.valid) & (hist.cbins >= tail_min)
    if not np.any(selected):
        return total / bins, total * 0
    tail = copy.copy(hist)
    tail.valid = selected
    chi2, ndf = objective.chi2_cost(tail, measured)
    if ndf < 1:
        raise FloatingPointError('No finite STAR tail residuals at the proposed amplitude point')
    return total / bins, chi2 / ndf


# Evaluate constrained pion fits with fixed native MC errors and autograd derivatives
class Refinement:
    # Keep the native reference and the common-event response at the same fitted point
    def __init__(self, responses, native, parameters, names, bounds, tail_min, tail_weight, output,
                 cms_weight=0.0, kaon_weight=0.0, tail_max=None):
        self.responses, self.native, self.parameters = responses, native, parameters
        self.names, self.bounds = names, np.asarray(bounds)
        self.tail_min, self.tail_weight, self.output = tail_min, tail_weight, output
        self.cms_weight, self.tail_max = cms_weight, tail_max
        self.kaon_weight = kaon_weight
        self.records, self.cached, self.value = [], None, None
        self.failures = []
        with torch.no_grad():
            self.anchor = [response.predict(parameters) for response in responses]
        self.baseline = None

    # Compute both the physics objective and the CMS agreement constraint
    def evaluate(self, unit):
        if self.cached is not None and np.array_equal(unit, self.cached):
            return self.value
        started = time.monotonic()
        point = self.bounds[:, 0] + np.asarray(unit) * (self.bounds[:, 1] - self.bounds[:, 0])
        (self.output / 'last_trial.json').write_text(json.dumps(
            dict(parameters=dict(zip(self.names, point.tolist(), strict=True))), indent=2) + '\n')
        x = torch.tensor(unit, dtype=torch.float64, requires_grad=True)
        costs, derivatives = [], []
        for index, response in enumerate(self.responses):
            values = x.new_tensor(self.bounds[:, 0]) + x * x.new_tensor(self.bounds[:, 1] - self.bounds[:, 0])
            parameters = self.parameters | dict(zip(self.names, values, strict=True))
            prediction = corrected_histograms(response.predict(parameters), self.anchor[index],
                                              self.native['mc'][index])
            total, tail = residuals(prediction, self.native['data'][index], self.tail_min)
            if not torch.isfinite(total) or not torch.isfinite(tail):
                raise FloatingPointError('Nonfinite histogram cost at the proposed amplitude point')
            total_gradient = (torch.autograd.grad(total, x, retain_graph=index == 1)[0]
                              if total.requires_grad else torch.zeros_like(x))
            tail_gradient = (torch.autograd.grad(tail, x)[0]
                             if index == 1 and tail.requires_grad else torch.zeros_like(x))
            costs.append((float(total.detach()), float(tail.detach())))
            derivatives.append((total_gradient.numpy(), tail_gradient.numpy()))
        cms, star, kaon, tail = costs[0][0], costs[1][0], costs[2][0], costs[1][1]
        objective = self.cms_weight * cms + star + self.kaon_weight * kaon + self.tail_weight * tail
        outputs = np.array((objective, cms, star, kaon, tail))
        gradients = np.stack((self.cms_weight * derivatives[0][0] + derivatives[1][0]
                              + self.kaon_weight * derivatives[2][0] + self.tail_weight * derivatives[1][1],
                              derivatives[0][0], derivatives[1][0], derivatives[2][0], derivatives[1][1]))
        if not np.all(np.isfinite(gradients)):
            raise FloatingPointError('Nonfinite amplitude gradient at the proposed point')
        self.value = outputs, gradients
        self.cached = np.array(unit, copy=True)
        record = dict(parameters=dict(zip(self.names, values.detach().tolist(), strict=True)),
                      objective=objective, cms=cms, star=star, kaon=kaon, tail=tail, elapsed=time.monotonic() - started)
        self.records.append(record)
        (self.output / 'evaluations.json').write_text(json.dumps(self.records, indent=2) + '\n')
        print(len(self.records), json.dumps(record), flush=True)
        return self.value

    # Reject nonfinite line search points without terminating the fit
    def safe_evaluate(self, unit):
        try:
            return self.evaluate(unit)
        except FloatingPointError as error:
            self.failures.append(dict(reason=str(error), unit=np.asarray(unit).tolist()))
            (self.output / 'failures.json').write_text(json.dumps(self.failures, indent=2) + '\n')
            self.cached = np.array(unit, copy=True)
            self.value = np.full(5, np.inf), np.zeros((5, len(unit)))
            print('Rejected amplitude point:', error, flush=True)
            return self.value

    # Compute the SLSQP objective and its gradient in normalized parameter coordinates
    def objective(self, unit):
        values, gradients = self.safe_evaluate(unit)
        return values[0], gradients[0]

    # Require CMS and both STAR comparisons to stay at least as good
    def constraints(self, unit):
        values, _ = self.safe_evaluate(unit)
        constraints = self.baseline[1:4] - values[1:4]
        return constraints if self.tail_max is None else np.append(constraints, self.tail_max - values[4])

    # Differentiate the three cross-section agreement constraints
    def jacobian(self, unit):
        _, gradients = self.safe_evaluate(unit)
        return -gradients[1:4] if self.tail_max is None else -gradients[1:5]


# Check directional derivatives on the actual screened bank before a global refinement
def check_gradient(fit, initial):
    direction = np.random.default_rng(29173).normal(size=len(initial))
    direction[(initial < 1e-4) | (initial > 1 - 1e-4)] = 0
    direction /= np.linalg.norm(direction)
    values, gradient = fit.evaluate(initial)
    step = 1e-6
    high = fit.evaluate(initial + step * direction)[0]
    low = fit.evaluate(initial - step * direction)[0]
    actual, expected = (high - low) / (2 * step), gradient @ direction
    error = np.abs(actual - expected) / np.maximum(1, np.abs(expected))
    report = dict(values=values.tolist(), finite_difference=actual.tolist(),
                  autograd=expected.tolist(), relative_error=error.tolist())
    (fit.output / 'gradient_check.json').write_text(json.dumps(report, indent=2) + '\n')
    if np.max(error) > 1e-4:
        raise ValueError('Screened refinement gradient fails the directional finite difference check')


# Rescale optimizer steps without changing any physical parameter bounds
def minimize_fit(fit, initial, scale, maxiter, method):
    # Map local optimizer coordinates onto the complete physical fit interval
    def coordinates(point):
        return np.clip(initial + scale * point, 0.0, 1.0)

    # Compute the scaled objective gradient for the same coherent predictions
    def objective(point):
        value, gradient = fit.objective(coordinates(point))
        return value, scale * gradient

    bounds = list(zip(-initial / scale, (1 - initial) / scale, strict=True))
    if method != 'slsqp':
        raise ValueError('Refinement requires SLSQP to enforce the data and tail constraints')
    result = minimize(objective, np.zeros_like(initial), jac=True, method='SLSQP', bounds=bounds,
                      constraints=[dict(type='ineq', fun=lambda point: fit.constraints(coordinates(point)),
                                        jac=lambda point: scale * fit.jacobian(coordinates(point)))],
                      options=dict(maxiter=maxiter, ftol=1e-8, disp=True))
    result.x = coordinates(result.x)
    return result


# Scale each coordinate by its local physical magnitude instead of its prior width
def relative_scale(names, bounds, initial):
    width = bounds[:, 1] - bounds[:, 0]
    physical = np.abs(bounds[:, 0] + initial * width)
    for index, name in enumerate(names):
        if name.endswith(('@PHASE', '@PROJECTIVE')):
            physical[index] = 1.0
    return np.clip(physical / width, 0.01, 1.0)


# Refine requested amplitude parameters while retaining the native reference comparisons
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tail-min', type=float, default=2.0)
    parser.add_argument('--tail-weight', type=float, default=1.0)
    parser.add_argument('--cms-weight', type=float, default=0.0)
    parser.add_argument('--kaon-weight', type=float, default=0.0)
    parser.add_argument('--method', choices=('slsqp',), default='slsqp')
    parser.add_argument('--tail-max', type=float, help='Optional upper constraint on the STAR tail chi2 per bin')
    parser.add_argument('--anchor', type=Path, help='Refinement point matching the supplied native histograms')
    parser.add_argument('--start-refinement', type=Path, help='Initial fit point with the native anchor kept fixed')
    parser.add_argument('--maxiter', type=int, default=60)
    parser.add_argument('--coordinate-scale', type=float, default=0.1)
    parser.add_argument('--relative-scale', action='store_true', help='Normalize steps to local parameter magnitudes')
    parser.add_argument('--parameters', choices=('pion', 'all'), default='pion')
    parser.add_argument('--check-gradient', action='store_true')
    parser.add_argument('--production-scale-max', type=float,
                        help='Optional upper bound for analytic hadronic FF_prod.Lambda2')
    args = parser.parse_args()
    if args.coordinate_scale <= 0:
        parser.error('--coordinate-scale must be positive')
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'arguments.json').write_text(json.dumps(vars(args), indent=2, default=str) + '\n')
    history = json.loads((args.run / 'history.json').read_text())
    best = min(history['trials'], key=lambda trial: trial['metrics']['chi2'])
    parameters = best['config']
    if args.anchor is not None:
        parameters = json.loads(args.anchor.read_text())['parameters']
    limits = {entry['name']: (entry['lower'], entry['upper']) for entry in history['parameter_space']}
    names = sorted(name for name in limits if args.parameters == 'all' or name.startswith('CON_GP|990:[211,211]:'))
    if args.production_scale_max is not None:
        for name in names:
            if name.endswith(':FF_prod.Lambda2') and not any(vector in name for vector in ('rho_770', 'phi_1020')):
                limits[name] = (limits[name][0], args.production_scale_max)
    bounds = np.asarray([limits[name] for name in names])
    controls = load_settings(campaign_source("tune-gpom-star-cms-ampfit") + "#/optimizer/ampfit")['bank']
    banks = load_banks(Path.cwd(), args.run, controls)
    responses = (banks if args.parameters == 'all' else
                 [PionResponse(bank, bank.initial | parameters) for bank in banks])
    parameters = banks[0].initial | parameters
    with args.native.open('rb') as stream:
        native = pickle.load(stream)['results']
    fit = Refinement(responses, native, parameters, names, bounds, args.tail_min, args.tail_weight, args.output,
                     cms_weight=args.cms_weight, kaon_weight=args.kaon_weight, tail_max=args.tail_max)
    initial = (np.array([parameters[name] for name in names]) - bounds[:, 0]) / np.diff(bounds)[:, 0]
    fit.baseline = fit.evaluate(initial)[0]
    if args.start_refinement is not None:
        start = json.loads(args.start_refinement.read_text())['parameters']
        initial = (np.array([start[name] for name in names]) - bounds[:, 0]) / np.diff(bounds)[:, 0]
    initial = np.clip(initial, 0, 1)
    if args.check_gradient:
        check_gradient(fit, initial)
    scale = args.coordinate_scale * (relative_scale(names, bounds, initial) if args.relative_scale else 1)
    result = minimize_fit(fit, initial, scale, args.maxiter, args.method)
    if not result.success or np.any(fit.constraints(result.x) < -1e-8):
        raise RuntimeError(f'Constrained refinement failed: {result.message}')
    updated = parameters | dict(zip(names, bounds[:, 0] + result.x * np.diff(bounds)[:, 0], strict=True))
    payload = dict(parameters=updated, active=names, baseline=fit.baseline.tolist(),
                   fitted=fit.evaluate(result.x)[0].tolist(), success=bool(result.success),
                   message=result.message, iterations=int(result.nit),
                   source_run=str(args.run), native=str(args.native),
                   tail_min=args.tail_min, tail_weight=args.tail_weight, cms_weight=args.cms_weight,
                   kaon_weight=args.kaon_weight, method=args.method, tail_max=args.tail_max, bounds=bounds.tolist())
    (args.output / 'refinement.json').write_text(json.dumps(payload, indent=2) + '\n')
    (args.output / 'points.json').write_text(json.dumps([updated], indent=2) + '\n')


if __name__ == '__main__':
    main()
