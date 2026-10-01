# Compare screened GP generator predictions with STAR and CMS data
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import copy
import json
import logging
import pickle
import shutil
from pathlib import Path

import numpy as np
import pyjson5
from core.stats import hist
from core.stats.objective import evaluate_cost_bundle
from core.tune import push
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.tunesetup import load_tunesetup

from submit import campaign_source


# Export the best evaluated point with the saved source tune and original card formatting
def export_tune(history, tune, root, directory):
    history = history.resolve()
    payload = json.loads(history.read_text())
    metric = payload['optimization']['cost']
    best = min(payload['trials'], key=lambda trial: trial['metrics'][metric])
    target = (root / 'modeldata' / tune).resolve()
    if target.parent != (root / 'modeldata').resolve() or target.exists():
        raise ValueError('--history requires a new tune name directly under modeldata')
    source = history.parent / 'results/amplitude/0/nominal/tune'
    rendered = {}
    for path in source.rglob('*.json'):
        relative = path.relative_to(source)
        template = root / 'modeldata/TUNE0' / relative
        original, prepared = template.read_text(), path.read_text()
        if original == prepared:
            continue
        rendered[target / relative], _ = push.render_json5_scalars(original, pyjson5.loads(prepared))
    shutil.copytree(source, target)
    push.atomic_write_texts(rendered)
    changes = GraniittiDriver().push_parameters(
        summary=best, target_path=str(target), cdir=str(root), options={}, confirm=lambda rows: True)
    (directory / 'parameters.json').write_text(json.dumps(best, indent=2) + '\n')
    (directory / 'tune_changes.json').write_text(json.dumps(changes, indent=2) + '\n')
    print('Exported', best['trial_id'], 'to', target, flush=True)


# Export the constrained continuum refinement through the normal parameter updater
def export_refinement(refinement, source_tune, tune, root, directory):
    payload = json.loads(refinement.read_text())
    target = (root / 'modeldata' / tune).resolve()
    source = (root / 'modeldata' / source_tune).resolve()
    if target.parent != (root / 'modeldata').resolve() or target.exists():
        raise ValueError('--refinement requires a new tune name directly under modeldata')
    if not source.is_dir():
        raise ValueError('--source-tune must identify the native reference tune')
    shutil.copytree(source, target)
    summary = {'config': {name: payload['parameters'][name] for name in payload['active']}}
    changes = GraniittiDriver().push_parameters(
        summary=summary, target_path=str(target), cdir=str(root), options={}, confirm=lambda rows: True)
    (directory / 'parameters.json').write_text(json.dumps(payload, indent=2) + '\n')
    (directory / 'tune_changes.json').write_text(json.dumps(changes, indent=2) + '\n')
    print('Exported refinement to', target, flush=True)


# Integrate generated histogram bins onto the fit bin edges
def fit_histograms(results, driver):
    mc = copy.deepcopy(results['mc'])
    for index, subsets in enumerate(driver.data):
        for subset, histograms in enumerate(subsets):
            for name, item in histograms.items():
                target = item['hdata']
                source = mc[index][subset][name]['hdata']
                groups = np.searchsorted(source.bins, target.bins)
                if np.any(groups >= len(source.bins)) or not np.allclose(source.bins[groups], target.bins):
                    raise ValueError('Fit bins must coincide with original histogram bin edges')
                bins, counts, errors = hist.rebin_histogram_groups(
                    source.bins, source.counts_scaled, source.errs_scaled, groups, True)
                mc[index][subset][name]['hdata'] = hist.hobj(
                    counts=counts, errs=errors, bins=bins, cbins=hist.edge2centerbins(bins),
                    valid=target.valid, binscale=1.0, density=False)
    return dict(mc=mc, data=driver.data, obs=driver.obs, datasets=driver.datasets)


# Save the standard histogram objective and a compact dataset comparison
def save_results(results, directory):
    directory.mkdir(parents=True, exist_ok=True)
    costs = evaluate_cost_bundle(results=results, selected_cost='chi2', cost_rho='quadratic',
                                 cost_avg='global-mean', covariance_payload=None, rngseed=29173, wasserstein_cache=None)
    with (directory / 'results.pkl').open('wb') as stream:
        pickle.dump(dict(results=results, costs=costs), stream)
    (directory / 'costs.json').write_text(json.dumps(costs, indent=2,
        default=lambda value: value.tolist() if isinstance(value, np.ndarray) else value.item()) + '\n')
    summary = []
    for index, dataset in enumerate(results['datasets']):
        name = dataset['sets'][0]['name']
        chi2 = sum(sum(row.values()) for row in costs['cost_arr']['chi2'][index])
        bins = sum(sum(row.values()) for row in costs['ndf_arr']['chi2'][index])
        summary.append(dict(dataset=name, chi2=chi2, bins=bins, chi2_per_bin=chi2 / bins))
    (directory / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(directory, json.dumps(summary), flush=True)
    return costs


# Generate independent events and retain histogram costs and comparison plots
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tune', default='TUNE0')
    parser.add_argument('--name', default='baseline')
    exports = parser.add_mutually_exclusive_group()
    exports.add_argument('--history', type=Path, help='Export the best history point to a new --tune before validation')
    exports.add_argument('--refinement', type=Path, help='Export a constrained refinement to a new --tune')
    parser.add_argument('--source-tune', help='Native reference tune required with --refinement')
    parser.add_argument('--output-dir', type=Path, default=Path('runs/icetune/tune-gpom-star-cms-ampfit/validation'))
    parser.add_argument('--reuse-results', type=Path, help='Evaluate fit bins from a saved original-bin generator comparison')
    parser.add_argument('--plots', action='store_true', help='Render the original-bin generator comparisons')
    parser.add_argument('--cms-events', type=int, default=300000)
    parser.add_argument('--star-events', type=int, default=30000)
    parser.add_argument('--max-time', type=int, default=14400, help='Generator time limit in seconds')
    args = parser.parse_args()
    if min(args.cms_events, args.star_events, args.max_time) <= 0:
        parser.error('Event counts and generator time limit must be positive')
    if args.refinement is not None and args.source_tune is None:
        parser.error('--refinement requires --source-tune')
    logging.basicConfig(level=logging.INFO)
    root = Path.cwd()
    directory = args.output_dir.resolve() / args.name
    directory.mkdir(parents=True, exist_ok=True)
    if args.history is not None:
        export_tune(args.history, args.tune, root, directory)
    if args.refinement is not None:
        export_refinement(args.refinement, args.source_tune, args.tune, root, directory)
    tunesetup = load_tunesetup(cdir=root, simdriver='graniitti',
                          name=campaign_source("tune-gpom-star-cms-ampfit"))
    for index, card in enumerate(tunesetup.datacards):
        card['nevents'] = args.cms_events if index == 0 else args.star_events
    original = copy.deepcopy(tunesetup.datacards)
    original[0]['datacard'] = 'icepack/SOFTCEP/CMS_2752118/dataset.json'
    driver = GraniittiDriver()
    driver.init_data(str(directory), original, 'default', str(root), pickle_dump=False)
    if args.reuse_results is None:
        mc = driver.compute(tunename=args.tune, datacards=original,
                            mc_steer={'tune_default': args.tune, 'tunesetup_name': args.name},
                            cdir=str(root), processes=3, max_t=args.max_time, rngseed=29173,
                            init_log_dir=str(directory / 'logs'))
        results = dict(mc=mc, data=driver.data, obs=driver.obs, datasets=driver.datasets)
    else:
        with args.reuse_results.open('rb') as stream:
            results = pickle.load(stream)['results']
    costs = save_results(results, directory)
    if args.plots:
        driver.render_trial_figures_to_dir(
            outputs=dict(results=results, valid_arr=costs['valid_arr'], likelihood=costs['likelihood'],
                         tunename=args.tune, trial_id=args.name),
            param=dict(run_name=args.name, cdir=str(root)), summary_payload={}, output_dir=str(directory / 'plots'))
    fit_cards = copy.deepcopy(tunesetup.datacards)
    fit_cards[0]['datacard'] = "icepack/SOFTCEP/CMS_2752118/dataset_icetune.json"
    fit_driver = GraniittiDriver()
    fit_driver.init_data(str(directory), fit_cards, 'default', str(root), pickle_dump=False)
    save_results(fit_histograms(results, fit_driver), directory.with_name(args.name + '_fitbins'))


if __name__ == '__main__':
    main()
