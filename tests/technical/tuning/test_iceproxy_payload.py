# Tests for iceproxy payloads and surrogate models
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
import json
import pickle
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from submit import CAMPAIGN_DIR
from submit import campaign as campaign_config
from submit.lxplus import outputs as run_outputs

ROOT = Path(__file__).resolve().parents[3]

from core.inference import gp as gp_model
from core.inference import mcmc
from core.stats import cov, hist
from core.stats import objective as icecost
from core.tune import core as icetune_main
from core.tune.backends import ray as icetune_ray  # noqa: E402
from core.tune.parameters import tools  # noqa: E402


# Load the extensionless iceproxy executable as a test module
def _load_iceproxy_module():
    from importlib import import_module

    return import_module('core.iceproxy')


# Load the extensionless icescape executable as a test module
def _load_icescape_module():
    from importlib import import_module

    return import_module('core.icescape')


# Select the replica backend through the common model option
def test_iceproxy_model_cli(monkeypatch, tmp_path):
    iceproxy = _load_iceproxy_module()
    run_directory = tmp_path / "runs" / "icetune" / "replica_cli"
    (run_directory / "results").mkdir(parents=True)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "iceproxy",
            "--input",
            str(run_directory),
            "--cached",
            "--mc-errors",
            "--posterior",
            "--profile",
            "true",
            "--plot-2d",
            "--profile-multistart",
            "--profile_workers",
            "5",
            "--profile_dpi",
            "280",
            "--no-gp-optimize",
            "--no-optimize",
            "--impact-mode",
            "fixed",
            "--plot-brand",
            "CUSTOM",
            "--no-predictive-plots",
        ],
    )

    args = iceproxy.parse_args()

    assert args.model == "gp"
    assert iceproxy.replica_backend_name(args.model) == "gp_pca"
    assert args.cached and args.mc_errors and args.posterior
    assert args.profile and args.plot_2d and args.profile_multistart
    assert args.profile_dpi == 280
    assert args.profile_workers == 5
    assert not args.gp_optimize and not args.optimize
    assert args.optimizer_backend == "torch"
    assert args.optimizer_lbfgs_maxiter == 200
    assert args.optimizer_restarts == 32
    assert args.optimizer_scan_points == 1024
    assert not hasattr(args, "impacts")
    assert args.impact_mode == "fixed"
    assert args.plot_brand == "CUSTOM"
    assert not args.predictive_plots
    assert args.input_mode == "full"
    assert args.run_name == "replica_cli"


# Build one minimal current-format full icetune pickle payload
def _write_trial(
    path: Path,
    x: float,
    data_counts: np.ndarray,
    data_errs: np.ndarray,
    *,
    data_covariance_mode: str = "diagonal",
    simdriver: str = "PANDORA",
) -> None:
    mc_counts = np.array([10.0 + 2.0 * x, 20.0 - x, 5.0 + x], dtype=float)
    mc_errs = np.array([0.2, 0.2, 0.1], dtype=float)
    data_rec = {
        "hdata": SimpleNamespace(counts_scaled=data_counts, errs_scaled=data_errs),
        "fitw": 1.0,
    }
    mc_rec = {
        "hdata": SimpleNamespace(counts_scaled=mc_counts, errs_scaled=mc_errs),
        "fitw": 1.0,
    }
    payload = {
        "replica_schema_version": 1,
        "trial_id": f"trial-{x:.2f}",
        "config": {"PARAM|x": x},
        "param": {
            "simdriver": simdriver,
            "plot_brand": None,
            "data_covariance_mode": data_covariance_mode,
        },
        "results": {
            "mc": [[{"obs": mc_rec}]],
            "data": [[{"obs": data_rec}]],
            "obs": None,
            "datasets": None,
        },
    }
    with open(path, "wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)


# Save one minimal data covariance beside full trial pickles
def _write_covariance(result_dir: Path, covariance: np.ndarray) -> None:
    size = len(covariance)
    data_stat_uncertainty = np.sqrt(np.diag(covariance))
    payload = {
        "schema_version": cov.DATA_COVARIANCE_SCHEMA,
        "method": "test",
        "asymmetric_rule": "mean_absolute_up_down",
        "layout": [
            {
                "dataset": 0,
                "subset": 0,
                "observable": "obs",
                "bins": list(range(size)),
                "start": 0,
                "stop": size,
            }
        ],
        "data_uncertainty_sources": [],
        "data_collective_sources": [],
        "mc_correlation_samples": [],
        "mc_correlation_scale": 1.0,
        "data_stat_uncertainty": data_stat_uncertainty,
        "data_count_stat_uncertainty": data_stat_uncertainty,
        **cov.sparse_covariance_arrays(covariance, "stat_local"),
        **cov.sparse_covariance_arrays(np.zeros_like(covariance), "syst_local"),
        **cov.sparse_covariance_arrays(np.zeros_like(covariance), "mc_cross"),
        **cov.covariance_decomposition(
            covariance,
            np.empty((size, 0), dtype=float),
        ),
    }
    cov.save_data_covariance(payload, output_dir=str(result_dir))


# Verify current full-pickle payloads flatten into replica tensors
def test_simu_pickle_payload(tmp_path):
    run_name = "replica_payload"
    result_dir = tmp_path / "runs" / "icetune" / run_name / "results"
    result_dir.mkdir(parents=True)

    data_counts = np.array([10.8, 19.6, 5.4], dtype=float)
    data_errs = np.array([1.0, 1.0, 0.5], dtype=float)
    for index, x in enumerate([0.0, 0.5, 1.0]):
        _write_trial(result_dir / f"TUNE_icetune_trial-{index:06d}.pkl", x, data_counts, data_errs)

    output = mcmc.collect_simu(run_name=run_name, cdir=str(tmp_path), max_trials=10)

    assert output["replica_input_schema_version"] == 1
    assert np.allclose(output["Y_data"], data_counts)
    assert np.allclose(output["E_data"], data_errs)
    assert output["simdriver"] == "PANDORA"
    assert output["plot_brand"] is None
    assert output["config_keys"] == ["PARAM|x"]
    assert output["X"].shape == (3, 1)
    assert output["Y"].shape == (3, 3)
    assert output["E"].shape == (3, 3)
    assert (result_dir / "combined.pkl").is_file()


# Reject changes in physical bin definitions across saved simulation trials
@pytest.mark.parametrize("change", ["edges", "valid", "fitw", "shape"])
def test_simu_changed_hist(tmp_path, change):
    run_name = "changed_bins"
    result_dir = tmp_path / "runs" / "icetune" / run_name / "results"
    result_dir.mkdir(parents=True)
    for index in range(2):
        path = result_dir / f"TUNE_icetune_{index}.pkl"
        _write_trial(path, float(index), np.ones(3), np.ones(3))
        with path.open("rb") as stream:
            payload = pickle.load(stream)
        for source in ("mc", "data"):
            record = payload["results"][source][0][0]["obs"]
            record["hdata"].bins = np.arange(4.0)
            if index == 1 and change == "edges":
                record["hdata"].bins *= 2.0
        if index == 1:
            record = payload["results"]["data"][0][0]["obs"]
            if change == "valid":
                record["hdata"].valid = np.array([True, False, True])
            elif change == "fitw":
                record["fitw"] = 2.0
            elif change == "shape":
                record["hdata"].counts_scaled = np.ones(2)
        with path.open("wb") as stream:
            pickle.dump(payload, stream)
    with pytest.raises(ValueError, match="histogram|aligned"):
        mcmc.collect_simu(run_name, str(tmp_path))


# Verify iceproxy Z uses the exact tune chi2 convention
def test_global_2nll_matches_tune_chi2():
    iceproxy = _load_iceproxy_module()
    y_data = np.array([10.0, 20.0], dtype=np.float64)
    e_data = np.array([1.0, 2.0], dtype=np.float64)
    y_hat = np.array([[11.0, 18.0]], dtype=np.float64)
    e_hat = np.array([[0.5, 1.0]], dtype=np.float64)

    var = e_data**2 + e_hat[0] ** 2
    expected = np.sum((y_data - y_hat[0]) ** 2 / var)
    z, z_err = iceproxy.global_2NLL(y_data, e_data, y_hat, e_hat)

    assert z.tolist() == pytest.approx([expected])
    assert z_err.shape == (1,)


# Verify the differentiable replica objective matches data covariance chi2
def test_correlated_torch_objective_matches_icetune():
    counts = torch.tensor([[11.0, 99.0, 18.0]], dtype=torch.float64, requires_grad=True)
    errors = torch.tensor([[0.5, 0.0, 1.0]], dtype=torch.float64)
    data = torch.tensor([10.0, 99.0, 20.0], dtype=torch.float64)
    data_errors = torch.tensor([1.0, 0.0, 2.0], dtype=torch.float64)
    weights = torch.tensor([1.5, 0.0, 0.7], dtype=torch.float64)
    mask = torch.tensor([True, False, True])
    indices = torch.tensor([0, 2], dtype=torch.long)
    covariance = torch.tensor([[1.0, 0.4], [0.4, 4.0]], dtype=torch.float64)

    objective = icecost.histogram_chi2_torch(
        counts=counts,
        errors=errors,
        data_counts=data,
        data_errors=data_errors,
        mask=mask,
        fit_weights=weights,
        data_total_covariance=covariance,
        covariance_indices=indices,
    )
    expected = icecost.correlated_chi2_arrays(
        mc_prediction=counts.detach().numpy()[0, [0, 2]],
        data_values=data.numpy()[[0, 2]],
        mc_stat_uncertainty=errors.numpy()[0, [0, 2]],
        fit_weights=weights.numpy()[[0, 2]],
        data_total_covariance=covariance.numpy(),
    )

    assert objective.detach().numpy().tolist() == pytest.approx([expected["chi2"]])
    objective.sum().backward()
    assert torch.all(torch.isfinite(counts.grad))


# Verify both full-histogram tools reconstruct the full-covariance likelihood
def test_full_hist_surfaces_data_cov(tmp_path):
    run_name = "correlated_full_surface"
    result_dir = tmp_path / "runs" / "icetune" / run_name / "results"
    result_dir.mkdir(parents=True)
    data_counts = np.array([10.8, 19.6, 5.4], dtype=float)
    data_errs = np.array([1.0, 1.0, 0.5], dtype=float)
    covariance = np.array(
        [
            [1.0, 0.4, 0.0],
            [0.4, 1.0, 0.1],
            [0.0, 0.1, 0.25],
        ]
    )
    for index, x in enumerate([0.0, 0.5, 1.0]):
        _write_trial(
            result_dir / f"TUNE_icetune_trial-{index:06d}.pkl",
            x,
            data_counts,
            data_errs,
            data_covariance_mode="full",
            simdriver="GRANIITTI",
        )
    _write_covariance(result_dir, covariance)
    args = SimpleNamespace(
        run_name=run_name,
        cdir=str(tmp_path),
        max_trials=10,
        cached=False,
    )

    replica_surface = _load_iceproxy_module().load_replica_surface(args)
    mesh_surface = _load_icescape_module().load_full_surface(args)

    assert replica_surface["covariance_mode"] == "full"
    assert mesh_surface["covariance_mode"] == "full"
    assert mesh_surface["Z"] == pytest.approx(replica_surface["Z"])
    assert mesh_surface["observable_chi2"].shape == (3, 1)
    assert mesh_surface["observable_chi2"][:, 0] == pytest.approx(mesh_surface["Z"])
    assert mesh_surface["observables"][0]["observable"] == "obs"


# Verify iceproxy identifies antipodal projective chart boundaries
def test_projective_feature_antipodes():
    iceproxy = _load_iceproxy_module()
    projective = [tools.projective_angle_key("RES|f0:GP:g_ls", index) for index in range(2)]
    names = [*projective, "REGGE|omega.GP"]
    bounds = np.array(
        [
            [-0.5 * np.pi, 0.5 * np.pi],
            [-0.5 * np.pi, 0.5 * np.pi],
            [0.0, 1.0],
        ],
        dtype=np.float64,
    )
    angle = 0.37
    values = np.array(
        [
            [angle, 0.5 * np.pi, 0.25],
            [-angle, -0.5 * np.pi, 0.25],
        ],
        dtype=np.float64,
    )
    topology = tools.build_parameter_topology(names)

    features = iceproxy.model_features(values, bounds, names, topology)

    assert features.shape == (2, 7)
    assert np.allclose(features[0], features[1], atol=1.0e-12)


# Verify the GP and RFF iceproxy executable paths on full trial pickles
@pytest.mark.parametrize("model", ["gp", "rff"])
def test_iceproxy_backend_basic(tmp_path, model):
    run_name = f"replica_{model}"
    backend_name = _load_iceproxy_module().replica_backend_name(model)
    result_dir = tmp_path / "runs" / "icetune" / run_name / "results"
    result_dir.mkdir(parents=True)

    data_counts = np.array([10.8, 19.6, 5.4], dtype=float)
    data_errs = np.array([1.0, 1.0, 0.5], dtype=float)
    for index, x in enumerate([0.0, 0.25, 0.5, 0.75, 1.0]):
        _write_trial(result_dir / f"TUNE_icetune_trial-{index:06d}.pkl", x, data_counts, data_errs)

    cmd = [
        sys.executable,
        "-m", "core.iceproxy",
        "--input",
        str(tmp_path / "runs" / "icetune" / run_name),
        "--model",
        model,
        "--profile",
        "false",
        "--quality_mode",
        "none",
        "--rff_features",
        "32",
        "--gp_opt_steps",
        "3",
        "--random_samples",
        "8",
        "--posterior",
        "--posterior_samples",
        "8",
        "--posterior_burnin",
        "2",
        "--impact-mode",
        "fixed",
    ]
    subprocess.run(cmd, check=True)

    output_runs = list((tmp_path / "figs" / "iceproxy").glob(f"{run_name}__*"))
    assert len(output_runs) == 1
    output_root = output_runs[0]
    backend_root = output_root / backend_name
    assert (output_root / "manifest.json").is_file()
    assert (backend_root / f"replica_{backend_name}_model.pt").is_file()
    assert (output_root / "summary.json").is_file()
    assert (backend_root / "parameters.json").is_file()
    for basis in ("optimizer", "physical"):
        assert (backend_root / basis / "parameter_uncertainties.png").is_file()
        assert (backend_root / basis / "parameter_uncertainties.pdf").is_file()
        assert (backend_root / basis / "impacts" / "fixed" / "manifest.json").is_file()
        assert (backend_root / basis / "posterior" / "posterior_samples.npz").is_file()
    posterior = backend_root / "optimizer" / "posterior"
    assert json.loads((posterior / "posterior_summary.json").read_text())["sampler"] == "reflective_hmc"
    with np.load(posterior / "posterior_diagnostics.npz") as diagnostics:
        assert len(diagnostics["accept_history"]) == 10
    summary = json.loads((output_root / "summary.json").read_text(encoding="utf-8"))
    parameters = json.loads((backend_root / "parameters.json").read_text(encoding="utf-8"))
    physical = json.loads((backend_root / "physical" / "parameters.json").read_text(encoding="utf-8"))
    np.testing.assert_allclose(physical["covariance"], parameters["covariance_matrix"])
    run_payload = json.loads((output_root / "manifest.json").read_text(encoding="utf-8"))
    assert run_payload["status"] == "completed"
    assert run_payload["input_run_name"] == run_name
    assert run_payload["arguments"]["impact_mode"] == "fixed"
    assert run_payload["arguments"]["optimizer_backend"] == "torch"
    assert parameters["optimization_diagnostics"]["backend"] == "torch"
    assert parameters["best_fit"]["uncertainty_method"] == "torch_autograd_hessian"
    assert (
        summary["objective_convention"]["definition"]
        == "fit-weighted sum of informative-bin chi2 terms"
    )
    assert summary["best_fit"]["schema_version"] == 1
    assert summary["best_fit"]["source"] == "surrogate"
    assert summary["best_fit"]["objective"]["name"] == "Z"
    assert set(summary["best_fit"]["parameters"]) == {"PARAM|x"}
    assert summary["best_fit"]["parameters"]["PARAM|x"]["uncertainty"] is not None


# Exercise the icescape Torch posterior after a torch fit on a real GP surrogate
def test_icescape_hmc_basic(tmp_path):
    run_dir = tmp_path / "runs" / "icetune" / "mesh_hmc"
    results = run_dir / "results"
    results.mkdir(parents=True)
    for index, x in enumerate(np.linspace(0.0, 1.0, 5)):
        _write_trial(results / f"TUNE_icetune_trial-{index:06d}.pkl", x,
                     np.array([10.8, 19.6, 5.4]), np.array([1.0, 1.0, 0.5]))
    subprocess.run([
        sys.executable, "-m", "core.icescape", "--input", str(run_dir), "--model", "gp",
        "--profile", "false", "--quality_mode", "none", "--gp_opt_steps", "3",
        "--optimizer_backend", "torch", "--posterior",
        "--posterior_samples", "8", "--posterior_burnin", "2", "--impact-mode", "fixed",
    ], check=True)
    output = next((tmp_path / "figs" / "icescape").glob("mesh_hmc__*/gp"))
    posterior = output / "optimizer" / "posterior"
    summary = json.loads((posterior / "posterior_summary.json").read_text())
    assert summary["sampler"] == "reflective_hmc"
    with np.load(posterior / "posterior_samples.npz") as samples:
        assert samples["samples"].shape == (8, 1)
        assert np.all((samples["samples_unit"] >= 0.0) & (samples["samples_unit"] <= 1.0))
        with np.load(output / "physical" / "posterior" / "posterior_samples.npz") as physical:
            np.testing.assert_allclose(physical["samples"], samples["samples"])


# Verify GP-PCA training uses validation-selected hyperparameters
def test_gp_pca_validation_selection():
    iceproxy = _load_iceproxy_module()
    x = np.linspace(0.0, 1.0, 8, dtype=np.float64)
    X_model = np.column_stack((x, 1.0 - x))
    Y_scaled = np.column_stack(
        (
            np.sin(x),
            np.cos(x),
            x,
            x**2,
        )
    )
    args = SimpleNamespace(
        validation_fraction=0.25,
        rngseed=1234,
        pca_variance=0.99,
        pca_max_modes=3,
        gp_scale=0.5,
        gp_noise=1e-4,
        gp_optimize=1,
        gp_opt_steps=3,
        gp_lr=0.05,
        gp_patience=3,
        gp_min_delta=0.0,
        device="cpu",
    )

    gp, pca, stats = iceproxy.train_gp_pca(X_model=X_model, Y_scaled=Y_scaled, args=args)
    assert gp.X_train.shape[0] == 8
    assert pca["n_modes"] >= 1
    assert len(stats["train_loss"]) == 3
    assert 0 < len(stats["validation_loss"]) <= len(stats["train_loss"])
    assert stats["best_validation_loss"] is not None
    assert stats["best_step"] is not None


# Verify an informative data bin cannot disappear from the replica objective
def test_fixed_mask_penalizes_zero_prediction():
    iceproxy = _load_iceproxy_module()
    z, _ = iceproxy.global_2NLL(
        Y_data=np.array([2.0]),
        E_data=np.array([0.5]),
        Y_hat=np.array([[0.0]]),
        E_hat=np.array([[0.1]]),
        mask=np.array([True]),
    )

    assert z[0] > 15.0


# Verify parameter covariance is propagated with exact Torch derivatives
def test_predictive_hist_error_autograd_jacobian():
    iceproxy = _load_iceproxy_module()

    # Compute an analytic histogram prediction with one fitted coordinate
    def predictor(values, include_model_variance=True):
        counts = torch.cat((2.0 * values, 3.0 * values), dim=1)
        errors = torch.full_like(counts, 0.1)
        model_variance = (
            torch.full_like(counts, 0.04) if include_model_variance else torch.zeros_like(counts)
        )
        return counts, errors, model_variance, torch.zeros_like(counts)

    prediction = iceproxy.predictive_histogram_summary(
        predict_histograms_torch=predictor,
        point=np.array([0.5]),
        covariance=np.array([[0.09]]),
        device="cpu",
        dtype=torch.float64,
    )

    assert prediction["parameter_sigma"] == pytest.approx([0.6, 0.9])
    assert prediction["model_count_sigma"] == pytest.approx([0.2, 0.2])
    assert np.all(prediction["total_sigma"] > prediction["parameter_sigma"])


# Verify the scalable GP shares one N by N factorization across outputs
def test_independent_gp_derivatives():
    X = torch.linspace(0.0, 1.0, 7, dtype=torch.float64)[:, None]
    Y = torch.cat((torch.sin(X), torch.cos(X), X**2), dim=1)
    gp = gp_model.GaussianProcess(
        X_train=X,
        Y_train=Y,
        device="cpu",
        dtype="float64",
        kernel="Matern",
        independent_outputs=True,
    )
    gp.freeze_for_inference()
    point = torch.tensor([[0.35]], dtype=torch.float64, requires_grad=True)
    mean, variance = gp.predict_marginals(point)
    (gradient,) = torch.autograd.grad(
        mean.sum() + variance.sum(),
        point,
        create_graph=True,
    )
    (hessian,) = torch.autograd.grad(gradient.sum(), point)

    assert gp.K.shape == (len(X), len(X))
    assert mean.shape == variance.shape == (1, Y.shape[1])
    assert torch.isfinite(gradient).all()
    assert torch.isfinite(hessian).all()


# Build one complete plotting record for predictive tune-style figures
def _plot_record(counts, errors, label, color):
    bins = np.array([0.0, 1.0, 2.0, 3.0])
    observable = {
        "xlabel": "$x$",
        "ylabel": "$d\\sigma/dx$",
        "units": {"x": "GeV", "y": "pb", "yden": "GeV"},
        "xlim": [0.0, 3.0],
        "ylim": None,
    }
    return {
        "hdata": hist.hobj(
            counts=np.asarray(counts, dtype=float),
            errs=np.asarray(errors, dtype=float),
            bins=bins,
            cbins=np.array([0.5, 1.5, 2.5]),
            binscale=1.0,
        ),
        "hfunc": "hist",
        "color": color,
        "label": label,
        "style": {"histtype": "step", "lw": 1.2},
        "obs": observable,
        "fitw": 1.0,
    }


# Verify lxplus publishes worker histograms and covariance that iceproxy can read
@pytest.mark.parametrize("covariance_mode", ["diagonal", "full"])
def test_lxplus_trial_outputs_load_iceproxy(tmp_path, covariance_mode):
    catalog = campaign_config.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml")
    environment = campaign_config.resolve(catalog, campaign_name="tune-gpom-res-con",)
    worker = tmp_path / "worker"
    head = tmp_path / "head"
    shared = tmp_path / "shared"
    run_name = "replica_transfer"
    fingerprint = "c" * 64
    head_run = head / "runs" / "icetune" / run_name
    result_dir = head_run / "results"
    result_dir.mkdir(parents=True)
    state = icetune_ray.BestState(experiment_dir=str(head_run), cdir=str(head), run_name=run_name, cost="chi2")
    param = {
        "cdir": str(worker), "run_name": run_name, "simdriver": "GRANIITTI",
        "pickle_dump": environment["FULL_OUTPUT"] == "1", "data_covariance_mode": covariance_mode,
    }
    measured = np.array([10.8, 19.6, 5.4])
    uncertainty = np.array([1.0, 1.0, 0.5])
    covariance = np.array([[1.0, 0.4, 0.0], [0.4, 1.0, 0.1], [0.0, 0.1, 0.25]])
    points = [0.0, 0.5, 1.0]
    expected = []
    for index, x in enumerate(points):
        counts = np.array([10.0 + 2.0 * x, 20.0 - x, 5.0 + x])
        errors = np.array([0.2, 0.2, 0.1])
        if covariance_mode == "full":
            cost = icecost.correlated_chi2_arrays(
                mc_prediction=counts, data_values=measured, mc_stat_uncertainty=errors,
                fit_weights=np.ones(3), data_total_covariance=covariance,
            )["chi2"]
        else:
            cost = icecost.global_histogram_chi2(measured, uncertainty, counts[None, :], errors[None, :])[0][0]
        expected.append(cost)
        outputs = {
            "trial_id": f"trial-{index:06d}", "config": {"PARAM|x": x}, "metrics": {"chi2": cost},
            "results": {
                "mc": [[{"obs": _plot_record(counts, errors, "MC", "red")}]],
                "data": [[{"obs": _plot_record(measured, uncertainty, "Data", "black")}]],
                "obs": None, "datasets": None,
            },
        }
        trial_dir = worker / str(index)
        relative = Path("results", icetune_main.trial_pickle_filename(outputs))
        descriptor = icetune_main.maybe_dump_trial_payload(outputs=outputs, param=param, destination=trial_dir / relative)
        assert descriptor is not None
        descriptor["relative_path"] = relative.as_posix()
        transfer = icetune_ray.collect_ray_trial_outputs(
            trial_dir=trial_dir, run_name=run_name, descriptor=descriptor, rendered_initial=False, rendered_best=False,
        )
        pending = state.publish_trial_outputs(
            trial_id=outputs["trial_id"], cost=cost, transfer=transfer, descriptor=descriptor,
            rendered_initial=False, rendered_best=False,
        )
        state.commit_trial_outputs(outputs["trial_id"], pending["output_token"])
        assert (result_dir / relative.name).read_bytes() == (trial_dir / relative).read_bytes()
    if covariance_mode == "full":
        _write_covariance(result_dir, covariance)
    payload = run_outputs.pack_runtime_outputs(figure_dir=None, result_files=sorted(result_dir.iterdir()), run_name=run_name)
    run_outputs.restore_outputs(payload=payload, shared_root=shared, run_name=run_name, campaign_fingerprint=fingerprint)
    campaign_run = run_outputs.campaign_run_path(shared, run_name, fingerprint)
    for source in result_dir.iterdir():
        assert (campaign_run / "results" / source.name).read_bytes() == source.read_bytes()
    args = SimpleNamespace(
        run_name=campaign_run.relative_to(shared / "runs/icetune").as_posix(),
        cdir=str(shared), max_trials=10, cached=False,
    )
    surface = _load_iceproxy_module().load_replica_surface(args)
    order = np.argsort(surface["X"][:, 0])
    np.testing.assert_allclose(surface["X"][order, 0], points)
    np.testing.assert_allclose(surface["Z"][order], expected)
    assert surface["covariance_mode"] == covariance_mode
    assert surface["Y_hat"].shape == (len(points), len(measured))


# Preserve already normalized predictions and errors for every template normalization
@pytest.mark.parametrize("density", [False, True])
@pytest.mark.parametrize("uncertainty", ["scaled", "shape"])
def test_predictive_hists_scaled_values(density, uncertainty):
    iceproxy = _load_iceproxy_module()
    original = hist.hobj(
        counts=np.array([4.0, 8.0]), errs=np.array([2.0, 3.0]),
        bins=np.array([0.0, 1.0, 3.0]), cbins=np.array([0.5, 2.0]),
        binscale=np.array([2.0, 0.0]), density=density, density_uncertainty=uncertainty,
        density_denominator_counts=np.array([8.0, 16.0]),
        density_denominator_errs=np.array([3.0, 4.0]),
    )
    template = {"mc": [[{"mass": {"hdata": original}}]]}
    manifest = [{"dataset": 0, "subset": 0, "observable": "mass", "start": 0, "stop": 2,
                 "valid": [True, True]}]
    prediction = {"counts": np.array([0.4, 0.3]), "total_sigma": np.array([0.07, 0.09])}
    result = iceproxy.fill_predictive_results(template, manifest, prediction)
    histogram = result["mc"][0][0]["mass"]["hdata"]
    np.testing.assert_allclose(histogram.counts_scaled, prediction["counts"])
    np.testing.assert_allclose(histogram.errs_scaled, prediction["total_sigma"])
    np.testing.assert_array_equal(original.counts, [4.0, 8.0])


# Verify predictive outputs reuse the tune plotting machinery
def test_predictive_hist_outputs(tmp_path):
    iceproxy = _load_iceproxy_module()
    template = {
        "mc": [[{"M": _plot_record([1.0, 1.0, 1.0], [0.1] * 3, "MC", "red")}]],
        "data": [[{"M": _plot_record([1.1, 1.8, 2.7], [0.2] * 3, "Data", "black")}]],
        "datasets": [{"type": "test", "sets": [{"name": "sample"}]}],
        "obs": None,
    }
    manifest = [
        {
            "dataset": 0,
            "subset": 0,
            "observable": "M",
            "start": 0,
            "stop": 3,
            "valid": [True, True, True],
            "fitw": 1.0,
            "ndf": 3,
        }
    ]
    prediction = {
        "counts": np.array([1.2, 1.9, 2.8]),
        "mc_errors": np.array([0.1, 0.1, 0.1]),
        "model_count_sigma": np.array([0.2, 0.2, 0.2]),
        "model_error_sigma": np.zeros(3),
        "parameter_sigma": np.array([0.3, 0.3, 0.3]),
        "total_sigma": np.sqrt(np.array([0.14, 0.14, 0.14])),
        "count_jacobian": np.ones((3, 1)),
    }
    surface = {
        "plot_results": template,
        "histogram_manifest": manifest,
        "plot_template_source": "trial.pkl",
        "param_names": ["x"],
    }
    args = SimpleNamespace(cdir=str(tmp_path), run_name="predictive")
    output_root = tmp_path / "figs" / "iceproxy" / "predictive" / "rff"
    result = iceproxy.render_predictive_histograms(
        surface=surface,
        prediction=prediction,
        point=np.array([0.5]),
        backend_name="pca_rff_ridge",
        output_root=output_root,
        args=args,
    )

    assert result["status"] == "rendered"
    assert Path(result["arrays"]).is_file()
    assert Path(result["summary"]).is_file()
    assert len(list(output_root.rglob("hplot__M.pdf"))) == 2


# Check covariance derivatives remain finite at repeated eigenvalues
def test_chi2_degenerate_eigenvalues():
    counts = torch.tensor([[1.0, 2.0]], dtype=torch.float64, requires_grad=True)
    errors = torch.ones_like(counts, requires_grad=True)

    # Compute the same correlated objective while varying predictions and MC uncertainty
    def objective(values, uncertainty):
        return icecost.histogram_chi2_torch(
            counts=values, errors=uncertainty, data_counts=torch.zeros(2, dtype=values.dtype),
            data_errors=torch.ones(2, dtype=values.dtype), mask=torch.ones(2, dtype=torch.bool),
            fit_weights=torch.ones(2, dtype=values.dtype), data_total_covariance=torch.eye(2, dtype=values.dtype),
            covariance_indices=torch.arange(2),
        )

    value = objective(counts, errors)
    gradient = torch.autograd.grad(value.sum(), errors)[0]
    torch.testing.assert_close(value, torch.tensor([2.5], dtype=counts.dtype))
    torch.testing.assert_close(gradient, torch.tensor([[-0.5, -2.0]], dtype=counts.dtype))
    assert torch.autograd.gradcheck(objective, (counts, errors))
    assert torch.autograd.gradgradcheck(objective, (counts, errors))


# Reject residuals outside the covariance support while retaining normalized shape residuals
def test_correlated_chi2_singular_support():
    counts = torch.tensor([[1.0, 1.0], [1.0, -1.0]], dtype=torch.float64)
    value = icecost.histogram_chi2_torch(
        counts=counts, errors=torch.zeros_like(counts), data_counts=torch.zeros(2, dtype=counts.dtype),
        data_errors=torch.ones(2, dtype=counts.dtype), mask=torch.ones(2, dtype=torch.bool),
        fit_weights=torch.ones(2, dtype=counts.dtype),
        data_total_covariance=torch.tensor([[1.0, -1.0], [-1.0, 1.0]], dtype=counts.dtype),
        covariance_indices=torch.arange(2),
    )
    assert torch.isinf(value[0])
    assert value[1].item() == pytest.approx(1.0)
