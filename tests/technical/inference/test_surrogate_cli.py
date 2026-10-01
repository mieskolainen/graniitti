# Technical tests for shared surrogate command-line runtime utilities
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from core.inference import surrogate_cli


# Compute the common quality-cut arguments used by surface filtering tests
def _quality_args(max_delta_z=2.0):
    return SimpleNamespace(
        quality_mode="delta",
        max_delta_Z=max_delta_z,
        quality_percentile=90.0,
        max_Z=float("inf"),
    )


# Compute the common optimizer arguments shared by icescape and iceproxy
def _optimizer_args(backend="torch"):
    return SimpleNamespace(
        optimizer_backend=backend,
        optimizer_lbfgs_maxiter=17,
        optimizer_restarts=2,
        optimizer_scan_points=8,
        rngseed=123,
        plot_brand="GRANIITTI",
    )


# Check input paths determine history and full-payload readers without mode flags
def test_history_or_full_run_input(tmp_path):
    history_run = tmp_path / "runs" / "icetune" / "history_run"
    history_path = history_run / "history.json"
    history_path.parent.mkdir(parents=True)
    history_path.write_text("{}", encoding="utf-8")
    history_args = surrogate_cli.normalize_input_args(
        SimpleNamespace(input=str(history_path), cdir=str(tmp_path)),
        replica=False,
    )

    assert history_args.input_mode == "history"
    assert history_args.history_path == str(history_path)
    assert history_args.run_name == "history_run"
    assert history_args.cdir == str(tmp_path)

    full_run = tmp_path / "runs" / "icetune" / "full_run"
    (full_run / "results").mkdir(parents=True)
    full_args = surrogate_cli.normalize_input_args(
        SimpleNamespace(input=str(full_run), cdir=str(tmp_path)),
        replica=True,
    )

    assert full_args.input_mode == "full"
    assert full_args.history_path is None
    assert full_args.run_name == "full_run"
    assert full_args.cdir == str(tmp_path)

    with pytest.raises(ValueError, match="full trial pickle"):
        surrogate_cli.normalize_input_args(
            SimpleNamespace(input=str(history_path), cdir=str(tmp_path)),
            replica=True,
        )


# Check invalid rows and quality rows are sliced identically across every tensor
def test_quality_filter_surface_row_alignment():
    surface = {
        "X": np.asarray([[0.0], [1.0], [np.nan], [3.0]]),
        "X_model": np.asarray([[10.0], [11.0], [12.0], [13.0]]),
        "Z": np.asarray([0.0, 1.0, 2.0, 10.0]),
        "aux": np.asarray([100, 101, 102, 103]),
    }
    filtered, threshold = surrogate_cli.quality_filter_surface(
        surface,
        finite_keys=("X", "X_model", "Z"),
        row_keys=("X", "X_model", "Z", "aux"),
        args=_quality_args(),
        minimum_rows=2,
    )
    np.testing.assert_array_equal(filtered["Z"], [0.0, 1.0])
    np.testing.assert_array_equal(filtered["X_model"].reshape(-1), [10.0, 11.0])
    np.testing.assert_array_equal(filtered["aux"], [100, 101])
    assert threshold == pytest.approx(2.0)


# Check the shared minimum-statistics guard is applied after all cuts
def test_quality_filter_surface_rejects_too_few_rows():
    surface = {
        "X": np.asarray([[0.0], [1.0], [2.0]]),
        "Z": np.asarray([0.0, 5.0, 10.0]),
    }
    with pytest.raises(RuntimeError, match="at least 3 trials"):
        surrogate_cli.quality_filter_surface(
            surface,
            finite_keys=("X", "Z"),
            row_keys=("X", "Z"),
            args=_quality_args(max_delta_z=0.1),
        )


# Persist surface dimensions through the real run manifest API
def test_record_surface_brand_and_dimensions(tmp_path):
    output = surrogate_cli.run_manifest.create_run_output(
        cdir=str(tmp_path), tool='icescape', tool_version=(0, 0, 1), input_run_name='surface', arguments={}, command=[])
    args = SimpleNamespace(output_root=output, plot_brand=None)
    surface = {'Z': np.array([1., 2., 3.]), 'param_names': ['LOOPSCREEN'],
               'simdriver': 'drivers/graniitti/driver.py'}
    surrogate_cli.record_surface(args, surface, inputs={'mode': 'history'}, histogram_bin_count=17)
    manifest = json.loads((output / 'manifest.json').read_text())
    assert manifest['inputs'] == {'mode': 'history'}
    assert manifest['surface'] == {'trial_count': 3, 'parameter_count': 1, 'histogram_bin_count': 17}


# Fit an analytic chi-square through the shared CLI and recover its Hessian covariance
@pytest.mark.parametrize('backend', ['torch'])
def test_optimize_surrogate_recovers_quadratic(tmp_path, backend):
    from core.stats.transform import ParameterTransform

    # Decode correlated sums and products through the shared uncertainty path
    def decode(values):
        return {"sum": values["x"] + values["y"], "product": values["x"] * values["y"]}

    # Evaluate scalar or batched points with either NumPy or Torch arrays
    def objective(values):
        return (values[..., 0] - .25)**2 + 4 * (values[..., 1] + .5)**2

    points = np.array([[-1., -1.], [0., 0.], [1., 1.]])
    result, best = surrogate_cli.optimize_surrogate(
        args=_optimizer_args(backend), X=points, Z=objective(points), param_names=['x', 'y'],
        func_hat=objective, bounds=np.array([[-1., 1.], [-1., 1.]]), output_root=tmp_path,
        realized_best=points[1], realized_Z=objective(points[1]), x0=np.array([-.5, .5]),
        objective_torch=objective, objective_batch_torch=objective, device='cpu', dtype=torch.float64,
        transform=ParameterTransform(['x', 'y'], decode, points[1]),
    )
    assert result.fmin.is_valid
    np.testing.assert_allclose(best, [.25, -.5], atol=1e-5)
    np.testing.assert_allclose(result.covariance, np.diag([1., .25]), atol=1e-5)
    payload = json.loads((tmp_path / 'parameters.json').read_text())
    physical = json.loads((tmp_path / 'physical' / 'parameters.json').read_text())
    jacobian = np.array([[1.0, 1.0], [best[1], best[0]]])
    np.testing.assert_allclose(physical['values'], [sum(best), np.prod(best)])
    np.testing.assert_allclose(physical['covariance'], jacobian @ np.asarray(result.covariance) @ jacobian.T)
    for name, expected in [('x', .25), ('y', -.5)]:
        assert payload['best_fit']['parameters'][name]['value'] == pytest.approx(expected, abs=1e-5)


# Render actual profile scans and honor the disabled profile path
@pytest.mark.parametrize('enabled', [False, True])
def test_render_profiles(tmp_path, enabled):
    args = SimpleNamespace(profile=enabled, ngrid=5, cmap='viridis_r', vmax=12., plot_2d=False,
                           profile_maxiter=20, profile_multistart=False, profile_workers=1,
                           profile_dpi=100, plot_brand='GRANIITTI')
    points = np.array([[0.], [.5], [1.]])
    rendered = surrogate_cli.render_profiles(
        args=args, X=points, Z=(points[:, 0] - .25)**2, param_names=['x'],
        func_hat=lambda values: (values[:, 0] - .25)**2, output_root=tmp_path,
        realized_best=points[0], surrogate_best=np.array([.25]), realized_Z=.25**2,
        surrogate_Z=0., bounds=np.array([[0., 1.]]),
    )
    assert rendered is enabled
    assert bool(list(tmp_path.rglob('hmesh__x.png'))) is enabled
