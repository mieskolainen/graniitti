# Resolve screened resonance and continuum contributions to the STAR pion mass tail
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import torch
from core.tune.drivers.graniitti.ampfit.amplitude import AmplitudeBank
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.tunesetup import load_tunesetup

from submit import campaign_source


# Integrate complete mass bins and their weighted MC uncertainties in each interval
def integrals(hist, edges):
    edges = np.asarray(edges, dtype=float)
    if len(edges) < 2 or not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise ValueError('Mass interval edges must be finite and strictly increasing')
    indices = np.abs(np.asarray(hist.bins)[:, None] - edges).argmin(axis=0)
    if not np.allclose(np.asarray(hist.bins)[indices], edges, rtol=1e-10, atol=1e-12):
        raise ValueError('Mass interval edges must coincide with histogram bin edges')
    rows = []
    for lower, upper in zip(edges[:-1], edges[1:], strict=True):
        selected = (hist.cbins >= lower) & (hist.cbins < upper)
        widths = np.diff(hist.bins)[selected]
        counts = np.asarray(hist.counts_scaled)[selected]
        errors = np.asarray(hist.errs_scaled)[selected]
        rows.append(dict(lower=lower, upper=upper, integral=float(counts @ widths),
                         error=float(np.linalg.norm(errors * widths))))
    return rows


# Contract a selected physical resonance or continuum before summing external helicities
def amplitude(bank, coefficients, selected):
    parts = []
    for start in range(0, len(bank.mass2), bank.batch_events):
        stop = start + bank.batch_events
        basis = bank.amplitudes[1:, start:stop][selected]
        parts.append((basis * coefficients[start:stop, selected].T[..., None]).sum(0))
    return torch.cat(parts)


# Save diagonal contributions and signed interference without dropping coherent phases
def inspect(bank, parameters, edges):
    with torch.no_grad():
        coefficients = bank.coefficients(bank.initial | parameters)
        size = len(bank.resonance.columns)
        groups = {name: [index for index, column in enumerate(bank.resonance.columns)
                         if column['resonance'] == name] for name in bank.resonance.resonances}
        groups['continuum'] = list(range(size, coefficients.shape[1]))
        amplitudes = {name: amplitude(bank, coefficients, selected) for name, selected in groups.items()}
        resonances = sum(value for name, value in amplitudes.items() if name != 'continuum')
        total = resonances + amplitudes['continuum']
        intensity = {name: value.abs().square().sum(1) for name, value in amplitudes.items()}
        intensity['resonances'] = resonances.abs().square().sum(1)
        intensity['total'] = total.abs().square().sum(1)
        result = {name: integrals(bank.weighted_histograms(value / bank.reference)[0]['M']['hdata'], edges)
                  for name, value in intensity.items()}
        # Difference coherent cross sections after histogramming to retain signed interference
        for name, positive, negative in (('resonance_continuum_interference', 'total', ['resonances', 'continuum']),
            ('resonance_interference', 'resonances', [key for key in groups if key != 'continuum']), ):
            result[name] = [{**row, 'error': None, 'integral': row['integral'] - sum(
                result[key][index]['integral'] for key in negative)} for index, row in enumerate(result[positive])]
        return result


# Inspect the same screened basis used by ampfit and compare with independent native histograms
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--native', required=True, type=Path)
    parser.add_argument('--parameters', type=Path, help='Refinement JSON instead of the source best history point')
    parser.add_argument('--mass-edges', required=True, type=float, nargs='+')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    history = json.loads((args.run / 'history.json').read_text())
    point = min(history['trials'], key=lambda trial: trial['metrics']['chi2'])['config']
    if args.parameters is not None:
        point = json.loads(args.parameters.read_text())['parameters']
    root, driver = Path.cwd(), GraniittiDriver()
    tunesetup = load_tunesetup(cdir=root, simdriver='graniitti',
                          name=campaign_source("tune-gpom-star-cms-ampfit"))
    driver.init_data(str(args.run), tunesetup.datacards, 'default', str(root), pickle_dump=False)
    controls = load_settings(campaign_source("tune-gpom-star-cms-ampfit") + "#/optimizer/ampfit")['bank']
    bank = AmplitudeBank(driver=driver, directory=args.run / 'results/amplitude/1/nominal',
                         obs=driver.obs[1], pid=driver.pid[1], cuts=driver.cuts[1], controls=controls)
    result = dict(bank=inspect(bank, point, args.mass_edges), parameters=point,
                  source_run=str(args.run), native_file=str(args.native))
    with args.native.open('rb') as stream:
        native = pickle.load(stream)['results']
    result['native'] = integrals(native['mc'][1][0]['M']['hdata'], args.mass_edges)
    result['data'] = integrals(native['data'][1][0]['M']['hdata'], args.mass_edges)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'parameters'}, indent=2))


if __name__ == '__main__':
    main()
