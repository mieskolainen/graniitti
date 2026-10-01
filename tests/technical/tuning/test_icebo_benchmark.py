# Tests for paired ICEBO optimizer benchmarks
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import numpy as np
import pytest
from core.tune.optimizers.icebo.benchmark import (
    ObjectiveEvaluation,
    _run_record,
    ackley,
    benchmark_optimizers,
    evaluate_observation,
    rastrigin,
    select_toy_problems,
    summarize_problems,
    toy_problems,
)


# Check the extended suite includes high-dimensional noisy objectives
def test_suite_contains_harder_noisy_problems():
    problems = {problem.name: problem for problem in toy_problems()}

    assert len(problems["ackley10_noisy"].bounds) == 10
    assert len(problems["rastrigin12_noisy"].bounds) == 12
    assert len(problems["levy16_noisy"].bounds) == 16
    assert all(
        problems[name].noise_std > 0.0
        for name in (
            "hartmann6_noisy",
            "ackley10_noisy",
            "rastrigin12_noisy",
            "levy16_noisy",
        )
    )


# Check new analytic objectives attain their declared origin minima
def test_new_objective_minima():
    assert ackley(np.zeros(10)) == pytest.approx(0.0, abs=1.0e-12)
    assert rastrigin(np.zeros(12)) == pytest.approx(0.0, abs=1.0e-12)


# Check noisy evaluations are paired, reproducible and heteroscedastic
def test_noise_stream_reproducible():
    problem = select_toy_problems(["ackley10_noisy"])[0]
    center = {f"x{index}": 0.0 for index in range(10)}
    corner = {f"x{index}": 5.0 for index in range(10)}

    first = evaluate_observation(problem, center, seed=17, evaluation_index=4)
    repeated = evaluate_observation(problem, center, seed=17, evaluation_index=4)
    changed = evaluate_observation(problem, center, seed=18, evaluation_index=4)
    outer = evaluate_observation(problem, corner, seed=17, evaluation_index=4)

    assert first == repeated
    assert first.observed != changed.observed
    assert outer.standard_deviation > first.standard_deviation


# Check deterministic objectives do not receive numerical pseudo-noise
def test_deterministic_observation_exact():
    problem = select_toy_problems(["branin2"])[0]
    config = {"x0": math.pi, "x1": 2.275}

    evaluation = evaluate_observation(problem, config, seed=17, evaluation_index=3)

    assert evaluation.standard_deviation == 0.0
    assert evaluation.observed == evaluation.latent


# Check noisy benchmark regret uses latent values rather than lucky fluctuations
def test_scores_latent_simple_regret():
    problem = select_toy_problems(["ackley10_noisy"])[0]
    evaluations = [
        ObjectiveEvaluation(latent=5.0, observed=-100.0, standard_deviation=1.0),
        ObjectiveEvaluation(latent=4.0, observed=-101.0, standard_deviation=1.0),
    ]

    record = _run_record("test", problem, 7, evaluations, 0.5, 1)

    assert record["best"] == 4.0
    assert record["best_observed"] == -101.0
    assert record["regret"] == 4.0
    assert record["normalized_regret"] == pytest.approx(0.8)


# Check per-problem summaries preserve optimizer wins and quantiles
def test_icebo_benchmark_problem_summary():
    runs = [
        {
            "problem": "p",
            "optimizer": optimizer,
            "seed": seed,
            "regret": regret,
            "normalized_regret": regret,
            "proposal_seconds": seconds,
        }
        for seed, values in ((1, (("icebo", 0.2, 1.0), ("hebo", 0.4, 2.0))),)
        for optimizer, regret, seconds in values
    ]

    summary = summarize_problems(runs)

    assert summary["p"]["icebo"]["wins"] == 1
    assert summary["p"]["hebo"]["wins"] == 0
    assert summary["p"]["icebo"]["median_normalized_regret"] == 0.2


# Run one actual paired noisy ICEBO and HEBO benchmark regression
@pytest.mark.parametrize("batch_size", [1, 2])
def test_runs_paired_noisy_comparison_against_hebo(batch_size):
    result = benchmark_optimizers(
        optimizers=["icebo", "hebo"],
        seeds=[11],
        budget=7,
        initial_points=4,
        fast=True,
        problem_names=["hartmann6_noisy"],
        batch_size=batch_size,
    )

    assert result["schema_version"] == 1
    assert {run["optimizer"] for run in result["runs"]} == {"icebo", "hebo"}
    assert all(run["initial_points"] == 4 for run in result["runs"])
    assert all(run["budget"] == len(run["latent_trace"]) == 7 for run in result["runs"])
    assert all(run["batch_size"] == batch_size for run in result["runs"])
    assert all(run["noise_std"] > 0.0 for run in result["runs"])
    assert all(math.isfinite(run["normalized_regret"]) for run in result["runs"])
    icebo, hebo = result["runs"]
    assert icebo["observed_trace"][:4] == hebo["observed_trace"][:4]
    assert set(hebo["numerics"]) == {
        "gp_fit_failures",
        "random_prediction_fallbacks",
    }


# Check parallel TPE includes pending trials and respects a partial final batch
def test_optuna_parallel_benchmark():
    pytest.importorskip("optuna")
    result = benchmark_optimizers(optimizers=["optuna"], seeds=[7], budget=8, initial_points=4,
                                  fast=True, problem_names=["hartmann6_noisy"], batch_size=3)
    run = result["runs"][0]
    assert run["budget"] == len(run["latent_trace"]) == 8
    assert np.isfinite(run["observed_trace"]).all()
    assert run["best"] <= run["latent_trace"][3]
