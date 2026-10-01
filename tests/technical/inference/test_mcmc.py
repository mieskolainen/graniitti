# Technical tests for shared MCMC numerical utilities
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
import pickle

import numpy as np
import pytest
import torch
from core.inference import mcmc


# Check cached trials obey the requested limit without slicing observed data
def test_collect_simu_limits_cached_trials(tmp_path):
    results = tmp_path / "runs" / "icetune" / "cached" / "results"
    results.mkdir(parents=True)
    payload = {
        "replica_input_schema_version": 1,
        **{name: np.arange(20).reshape(10, 2) for name in ("X", "Y", "E")},
        "Y_data": np.arange(5), "E_data": np.ones(5), "config_keys": ["a", "b"],
    }
    with (results / "combined.pkl").open("wb") as stream:
        pickle.dump(payload, stream)
    output = mcmc.collect_simu("cached", str(tmp_path), max_trials=2, use_cached=True)
    for name in ("X", "Y", "E"):
        np.testing.assert_array_equal(output[name], payload[name][:2])
    np.testing.assert_array_equal(output["Y_data"], payload["Y_data"])
    larger = mcmc.collect_simu("cached", str(tmp_path), max_trials=10, use_cached=True)
    assert larger["X"].shape[0] == 10
    with pytest.raises(ValueError, match="positive"):
        mcmc.collect_simu("cached", str(tmp_path), max_trials=0, use_cached=True)


# Integrate one scalar reflective trajectory through explicit collision times
def _scalar_reflection(position, momentum, duration, lower, upper):
    remaining = float(duration)
    position = float(position)
    momentum = float(momentum)
    while remaining > 1.0e-12 and momentum != 0.0:
        boundary = upper if momentum > 0.0 else lower
        hit_time = (boundary - position) / momentum
        if hit_time >= remaining:
            position += remaining * momentum
            break
        position = boundary
        remaining -= max(hit_time, 0.0)
        momentum = -momentum
    return position, momentum


# Check vectorized reflection against explicit trajectories with many collisions
def test_reflection_collision_integrator():
    rng = np.random.default_rng(104729)
    lower = rng.uniform(-3.0, -0.2, size=128)
    upper = lower + rng.uniform(0.1, 3.0, size=128)
    position = rng.uniform(lower, upper)
    momentum = rng.uniform(-5.0, 5.0, size=128)
    duration = 4.25

    reflected, reflected_momentum = mcmc.reflective_position_update(
        torch.tensor(position),
        torch.tensor(momentum),
        duration,
        torch.tensor(lower),
        torch.tensor(upper),
    )
    expected = [
        _scalar_reflection(x, p, duration, low, high)
        for x, p, low, high in zip(position, momentum, lower, upper, strict=True)
    ]
    np.testing.assert_allclose(reflected.numpy(), [item[0] for item in expected], atol=1.0e-12)
    np.testing.assert_allclose(
        reflected_momentum.numpy(),
        [item[1] for item in expected],
        atol=1.0e-12,
    )


# Check exact boundary handling and invalid boundary rejection
def test_reflection_boundaries():
    position = torch.tensor([0.5, 0.5, 0.5, 0.5], dtype=torch.float64)
    momentum = torch.tensor([1.0, -1.0, 1.0, -1.0], dtype=torch.float64)
    reflected, reflected_momentum = mcmc.reflective_position_update(
        position,
        momentum,
        1.5,
        torch.zeros(4, dtype=torch.float64),
        torch.ones(4, dtype=torch.float64),
    )
    torch.testing.assert_close(
        reflected,
        torch.tensor([0.0, 1.0, 0.0, 1.0], dtype=torch.float64),
    )
    torch.testing.assert_close(
        reflected_momentum,
        torch.tensor([-1.0, 1.0, -1.0, 1.0], dtype=torch.float64),
    )
    rounded_position, rounded_momentum = mcmc.reflective_position_update(
        torch.tensor([0.1], dtype=torch.float64),
        torch.tensor([0.2], dtype=torch.float64),
        1.0,
        torch.zeros(1, dtype=torch.float64),
        torch.tensor([0.3], dtype=torch.float64),
    )
    assert torch.equal(rounded_position, torch.tensor([0.3], dtype=torch.float64))
    assert torch.equal(rounded_momentum, torch.tensor([0.2], dtype=torch.float64))
    with pytest.raises(ValueError, match="lower < upper"):
        mcmc.reflective_position_update(
            position[:1],
            momentum[:1],
            1.0,
            torch.ones(1),
            torch.zeros(1),
        )


# Check vectorized posterior quantiles against NumPy definitions
def test_analyze_samples_quantiles_and_cov():
    samples = torch.tensor(
        [[0.0, 3.0], [1.0, 2.0], [2.0, 1.0], [3.0, 0.0]],
        dtype=torch.float64,
    )
    result = mcmc.analyze_samples(samples, config_keys=["a", "b"])
    np.testing.assert_allclose(result["median"], np.percentile(samples.numpy(), 50.0, axis=0))
    np.testing.assert_allclose(
        result["ci68"],
        np.percentile(samples.numpy(), [16.0, 84.0], axis=0).T,
    )
    np.testing.assert_allclose(result["covariance"], np.cov(samples.numpy(), rowvar=False))


# Check bounded HMC samples the Gaussian measure with the correct mean and covariance
def test_hmc_sampler_bounded_gaussian():
    torch.manual_seed(271828)
    warmup, retained = 512, 4096
    total = warmup + retained

    # Compute a differentiable standard-normal log posterior
    def log_posterior(position):
        return -0.5 * torch.sum(position.square()), True

    samples, diagnostics = mcmc.hmc_sampler(
        initial_position=torch.tensor([0.25, -0.5], dtype=torch.float64),
        n_samples=total,
        step_size=0.2,
        n_leapfrog_steps=5,
        log_posterior_fn=log_posterior,
        adapt_step_size=False,
        lower_bounds=torch.full((2,), -2.0, dtype=torch.float64),
        upper_bounds=torch.full((2,), 2.0, dtype=torch.float64),
        device="cpu",
    )
    assert samples.shape == (total, 2)
    assert torch.all(torch.isfinite(samples))
    assert torch.all((samples >= -2.0) & (samples <= 2.0))
    assert len(diagnostics["accept_history"]) == total
    assert len(diagnostics["energy_history"]) == total
    assert 0.0 < diagnostics["acceptance_rate"] <= 1.0

    # Integrate x^2 exp(-x^2/2) over [-2, 2] using its boundary term
    variance = 1.0 - 4.0 * math.exp(-2.0) / (math.sqrt(2.0 * math.pi) * math.erf(math.sqrt(2.0)))
    x = samples[warmup:].numpy()
    moments = np.column_stack((x, x**2, x[:, 0] * x[:, 1]))
    # Blocks longer than the HMC correlation time retain the Markov chain uncertainty
    batches = moments.reshape(64, 64, 5).mean(axis=1)
    error = batches.std(axis=0, ddof=1) / math.sqrt(len(batches))
    expected = np.array([0.0, 0.0, variance, variance, 0.0])
    assert np.all(error > 0.0)
    assert np.all(np.abs(batches.mean(axis=0) - expected) <= 5.0 * error)


# Compare one Gaussian HMC trajectory and Metropolis decision with the closed oscillator map
@pytest.mark.parametrize('step', [0.12, 0.6, 1.9])
@pytest.mark.parametrize('leaps', [1, 4])
def test_hmc_gaussian_hamiltonian(step, leaps):
    initial = torch.tensor([0.25, -0.5], dtype=torch.float64)
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(31)
        momentum = torch.randn_like(initial).numpy()
        uniform = torch.rand((), dtype=initial.dtype).item()
        theta = 2.0 * np.arcsin(step / 2.0)
        omega = np.sqrt(1.0 - step**2 / 4.0)
        cosine, sine = np.cos(leaps * theta), np.sin(leaps * theta)
        position = cosine * initial.numpy() + sine * momentum / omega
        final_momentum = -omega * sine * initial.numpy() + cosine * momentum
        energy = 0.5 * (position @ position + final_momentum @ final_momentum)
        delta = energy - 0.5 * (initial.numpy() @ initial.numpy() + momentum @ momentum)
        accepted = np.log(uniform) < min(0.0, -delta)
        torch.manual_seed(31)
        samples, diagnostics = mcmc.hmc_sampler(
            initial_position=initial, n_samples=1, step_size=step, n_leapfrog_steps=leaps,
            log_posterior_fn=lambda point: (-0.5 * point.square().sum(), True), adapt_step_size=False,
        )
    np.testing.assert_allclose(samples[0].numpy(), position if accepted else initial.numpy(), atol=1e-12)
    assert diagnostics['energy_history'] == pytest.approx([energy], rel=1e-12)
    assert diagnostics['accept_history'] == [float(accepted)]


# Check out-of-support HMC proposals count as rejections during adaptation
def test_hmc_invalid_trajectories_reduce_step_size():
    torch.manual_seed(271829)

    # Restrict the posterior to a narrow interval around the initial point
    def log_posterior(position):
        return -0.5 * position.square().sum(), bool(torch.all(position.abs() < 1e-12))

    samples, diagnostics = mcmc.hmc_sampler(
        initial_position=torch.zeros(1, dtype=torch.float64), n_samples=5,
        step_size=1.0, n_leapfrog_steps=2, log_posterior_fn=log_posterior,
    )
    assert torch.allclose(samples, torch.zeros_like(samples))
    assert diagnostics["accept_history"] == [0.0] * 5
    assert len(diagnostics["energy_history"]) == 5
    assert diagnostics["final_step_size"] < 1.0
    assert all(diagnostics["divergent_history"])
    with pytest.raises(ValueError, match="invalid log posterior"):
        mcmc.hmc_sampler(
            initial_position=torch.ones(1, dtype=torch.float64), n_samples=5,
            step_size=1.0, n_leapfrog_steps=2, log_posterior_fn=log_posterior,
        )


# Check short chains never request unavailable autocorrelation lags
@pytest.mark.parametrize("size", [1, 3, 12])
def test_autocorrelation_short_chain(size):
    values = np.arange(size, dtype=float)
    lags, correlation = mcmc.compute_autocorrelation(values)
    np.testing.assert_array_equal(lags, np.arange(size))
    assert correlation.shape == (size,)
    if size > 1:
        centered = values - values.mean()
        expected = [
            np.dot(centered[:size - lag], centered[lag:]) / (size - lag) / np.var(values)
            for lag in lags
        ]
        np.testing.assert_allclose(correlation, expected)
    else:
        assert np.isnan(correlation[0])


# Match direct lag products for a long chain without quadratic full correlation work
def test_long_chain_autocorrelation():
    values = np.random.default_rng(29).normal(size=32768)
    values[1:] += 0.5 * values[:-1]
    lags, actual = mcmc.compute_autocorrelation(values, max_lag=31)
    centered = values - values.mean()
    expected = np.array([
        np.dot(centered[:len(values) - lag], centered[lag:]) / (len(values) - lag)
        for lag in lags
    ])
    np.testing.assert_allclose(actual, expected / expected[0], atol=1.0e-14)


# Reject impossible initial box states before integrating a Hamiltonian trajectory
@pytest.mark.parametrize("lower, upper, start", [([0.0], [1.0], [2.0]), ([0.0], [math.inf], [0.5]),
    ([0.0], [0.0], [0.0]), ([0.0, 0.0], [1.0, 1.0], [0.5])])
def test_hmc_invalid_box_init(lower, upper, start):
    # Supply a valid Gaussian posterior independent of boundary validation
    def posterior(position):
        return -0.5 * position.square().sum(), True

    with pytest.raises(ValueError):
        mcmc.hmc_sampler(
            torch.tensor(start), n_samples=1, step_size=0.1, n_leapfrog_steps=1,
            log_posterior_fn=posterior, lower_bounds=torch.tensor(lower), upper_bounds=torch.tensor(upper),
        )
