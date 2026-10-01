# Test for shared posterior sampling
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pytest
import torch
from core.inference import mcmc


# Check saved physical coordinates and the Gaussian log posterior agree with the sampled points
@pytest.mark.parametrize('burnin', [0, 5])
def test_run_posterior_sampling_writes_samples(tmp_path, burnin):
    # Evaluate the same Gaussian objective on Torch and NumPy coordinates
    def objective(X):
        return ((X[..., 0] - 0.25) / 0.2) ** 2 + ((X[..., 1] + 0.5) / 0.3) ** 2

    output_dir = tmp_path / "posterior"
    samples, summary = mcmc.run_posterior_sampling(
        objective_torch=objective,
        bounds=np.array([[-1.0, 1.0], [-1.0, 1.0]], dtype=np.float64),
        param_names=["PARAM|x", "PARAM|y"],
        output_dir=output_dir,
        start=np.array([0.25, -0.5], dtype=np.float64),
        rngseed=123,
        n_samples=20,
        burnin=burnin,
        thin=3,
        step=0.05,
    )

    assert samples.shape == (20, 2)
    assert summary["n_samples"] == 20
    assert summary["sampler"] == "reflective_hmc"
    assert np.all(np.abs(samples) <= 1.0)
    with np.load(output_dir / "posterior_diagnostics.npz") as diagnostics:
        assert len(diagnostics["accept_history"]) == burnin + 60
        np.testing.assert_allclose(diagnostics["log_posterior_history"][burnin::3], -0.5 * objective(samples))
    with np.load(output_dir / "posterior_samples.npz") as saved:
        np.testing.assert_array_equal(saved["samples"], samples)
        np.testing.assert_allclose(saved["samples_unit"], (samples + 1.0) / 2.0)
        np.testing.assert_allclose(saved["log_posterior"], -0.5 * objective(samples))
    assert (output_dir / "posterior_summary.json").exists()
    assert (output_dir / "posterior__PARAM_x.png").exists()


# A zero integration step would falsely report perfect acceptance without sampling
@pytest.mark.parametrize("step", [0.0, -0.1, np.nan, np.inf])
def test_posterior_rejects_invalid_step(tmp_path, step):
    # Evaluate a one-dimensional Gaussian objective
    def objective(values):
        return values[..., 0] ** 2

    with pytest.raises(ValueError, match="HMC step"):
        mcmc.run_posterior_sampling(
            objective, np.array([[-1.0, 1.0]]), ["x"], tmp_path, np.array([0.0]), n_samples=1, step=step,
        )


# Check warmup freezes the averaged step and local generators preserve the global RNG
def test_hmc_warmup_and_local_rng():
    state = torch.random.get_rng_state().clone()
    options = dict(
        initial_position=torch.zeros(2, dtype=torch.float64), n_samples=70, step_size=0.1,
        n_leapfrog_steps=4, adaptation_window=50,
        log_posterior_fn=lambda q: (-0.5 * q.square().sum(), True),
    )
    first, diagnostics = mcmc.hmc_sampler(**options, generator=torch.Generator().manual_seed(137))
    second, _ = mcmc.hmc_sampler(**options, generator=torch.Generator().manual_seed(137))
    assert torch.equal(torch.random.get_rng_state(), state)
    torch.testing.assert_close(first, second, rtol=0, atol=0)
    np.testing.assert_allclose(diagnostics['step_size_history'][50:], diagnostics['final_step_size'])
    assert 0.0 < np.mean(diagnostics['accept_prob_history'][50:]) < 1.0


# Verify leapfrog reversibility, volume preservation and second order convergence on a correlated density
@pytest.mark.parametrize('bounded', [False, True])
def test_leapfrog_geometry(bounded):
    precision = torch.tensor([[2.0, 0.6], [0.6, 1.0]], dtype=torch.float64)
    initial = torch.tensor([0.2, -0.3, 0.8, -0.5], dtype=torch.float64)
    lower, upper = (torch.full((2,), -0.4), torch.full((2,), 0.4)) if bounded else (None, None)

    # Compute the correlated Gaussian target using the real autograd integrator
    def logp(q):
        return -0.5 * q @ precision @ q, True

    # Integrate a fixed duration from phase space coordinates
    def trajectory(z, steps=8):
        q, p, _, _ = mcmc.leapfrog(z[:2], z[2:], -precision @ z[:2], 0.8 / steps, steps, logp, lower, upper)
        return torch.cat((q, p))

    final = trajectory(initial)
    reverse = trajectory(torch.cat((final[:2], -final[2:])))
    torch.testing.assert_close(reverse, torch.cat((initial[:2], -initial[2:])), rtol=1e-12, atol=1e-12)
    displacement = torch.eye(4, dtype=torch.float64) * 1e-6
    jacobian = torch.stack([(trajectory(initial + d) - trajectory(initial - d)) / 2e-6 for d in displacement], dim=1)
    assert torch.linalg.det(jacobian).abs().item() == pytest.approx(1.0, abs=1e-8)
    if not bounded:
        generator = torch.zeros((4, 4), dtype=torch.float64)
        generator[:2, 2:] = torch.eye(2)
        generator[2:, :2] = -precision
        exact = torch.linalg.matrix_exp(0.8 * generator) @ initial
        errors = [torch.linalg.vector_norm(trajectory(initial, steps) - exact).item() for steps in [8, 16]]
        assert errors[0] / errors[1] == pytest.approx(4.0, rel=0.02)


# Recover correlated Gaussian moments after adaptation with Markov chain uncertainty from batches
def test_hmc_correlated_gaussian():
    covariance = torch.tensor([[1.0, 0.8], [0.8, 1.0]], dtype=torch.float64)
    precision = torch.linalg.inv(covariance)
    chain, diagnostics = mcmc.hmc_sampler(
        torch.zeros(2, dtype=torch.float64), n_samples=256 + 4096, step_size=0.1,
        n_leapfrog_steps=7, log_posterior_fn=lambda q: (-0.5 * q @ precision @ q, True),
        adaptation_window=256, step_jitter=0.2, generator=torch.Generator().manual_seed(139),
    )
    q = chain[256:].numpy()
    moments = np.column_stack((q, q**2, q[:, 0] * q[:, 1]))
    batches = moments.reshape(64, 64, 5).mean(axis=1)
    errors = batches.std(axis=0, ddof=1) / np.sqrt(len(batches))
    assert np.all(np.abs(batches.mean(axis=0) - [0, 0, 1, 1, 0.8]) < 5 * errors)
    assert not any(diagnostics['divergent_history'][256:])


# Reject non-finite gradients and divergent energy errors without losing the accepted state
def test_hmc_numerical_failures():
    for logp in [lambda q: (-torch.sqrt(q.abs()).sum(), True), lambda q: (-q.square().sum(), True)]:
        chain, diagnostics = mcmc.hmc_sampler(
            torch.ones(1, dtype=torch.float64), n_samples=4, step_size=1e5,
            n_leapfrog_steps=2, log_posterior_fn=logp, adapt_step_size=False,
            generator=torch.Generator().manual_seed(140),
        )
        assert torch.isfinite(chain).all()
        assert all(diagnostics['divergent_history'])
        torch.testing.assert_close(chain, torch.ones_like(chain))
    with pytest.raises(ValueError, match='gradient'):
        mcmc.hmc_sampler(torch.zeros(1), 1, 0.1, 1, lambda q: (-torch.sqrt(q.abs()).sum(), True))


# Count a finite-density proposal with an infinite autograd gradient as a rejection
def test_hmc_rejects_singular_gradient():
    momentum = torch.randn((), dtype=torch.float64, generator=torch.Generator().manual_seed(143)).item()
    samples, diagnostics = mcmc.hmc_sampler(
        torch.zeros(1, dtype=torch.float64), 1, 1.0 / abs(momentum), 1,
        lambda q: (torch.sqrt(1.0 - q.square()).sum(), True), adapt_step_size=False,
        generator=torch.Generator().manual_seed(143),
    )
    torch.testing.assert_close(samples, torch.zeros_like(samples))
    assert diagnostics['divergent_history'] == [True]
    assert diagnostics['accept_prob_history'] == [0.0]


# Keep adaptation finite on a uniform target whose Hamiltonian trajectories have unit acceptance
def test_hmc_uniform_box():
    samples, diagnostics = mcmc.hmc_sampler(
        torch.zeros(1, dtype=torch.float64), 256 + 1024, 0.1, 3, lambda q: (q.sum() * 0.0, True),
        adaptation_window=256, step_size_bounds=(1e-10, 1.0),
        lower_bounds=torch.zeros(1), upper_bounds=torch.ones(1), generator=torch.Generator().manual_seed(145),
    )
    q = samples[256:, 0].numpy()
    assert np.all((q >= 0.0) & (q <= 1.0))
    assert diagnostics['final_step_size'] <= 1.0
    assert diagnostics['acceptance_rate'] == pytest.approx(1.0)
    assert np.mean(q) == pytest.approx(0.5, abs=5 * np.sqrt(1.0 / (12 * len(q))))
    assert np.mean(q**2) == pytest.approx(1.0 / 3.0, abs=5 * np.sqrt(4.0 / (45 * len(q))))
