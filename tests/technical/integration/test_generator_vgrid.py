# Unit test MC integration & grids
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
#
# Run with: pytest tests/technical/integration/test_generator_vgrid.py --run-integration -s -rP

import json
import math

import json5
import matplotlib.pyplot as plt
import numpy as np
import pytest

from tests.technical.support.commands import execute
from tests.technical.support.output import PhysicsOutput

OUTPUT = PhysicsOutput("generator_vgrid")


# Compare multithreaded phase-space integration across mappings and final states
@pytest.mark.parametrize("FS", ["pi+ pi-", "K+ K-", "pi+ pi- pi+ pi-"])
def test_multicore(FS: str, physics_screening):
    """
    Test multithreaded sampling consistency of different amplitudes.
    Should yield the same results, if no multithreading related problems.
    """
    print(__name__ + ": Generating events ...")

    sigma = {"F": [], "C": []}
    errors = {"F": [], "C": []}

    # Let the unchanged 1% precision target determine the integration statistics
    for CORES in [1, 4]:
        for index, PS in enumerate(sigma):
            cmd = (
                f'./bin/gr -i ./tests/physics/studies/phase_space/gencard_four_body.json '
                f'-p "GP[CON]<{PS}> -> {FS}" -h 0 -n 0 -c {CORES} '
                f'-o multicore_{PS}_{CORES} -g VEGAS -r {27183 + 7919 * CORES + 104729 * index} '
                '--set INTEGRATOR.min_samples=10000'
            )
            physics_screening.execute(cmd)

            ## Read integration results
            with open(f"./vgrid/multicore_{PS}_{CORES}.vgrid") as f:
                data = json.load(f)
                sigma[PS].append(data["STAT"]["sigma"])
                errors[PS].append(data["STAT"]["sigma_err"])
                assert 0.0 < data["STAT"]["sigma_err"] <= 0.01 * sigma[PS][-1]

    # To Numpy
    for PS in sigma:
        sigma[PS] = np.array(sigma[PS])

    reference = sigma["F"][0]
    assert np.isfinite(reference)
    assert reference > 0.0
    for mode, values in sigma.items():
        assert np.isfinite(values).all()
        for rate, error in zip(values, errors[mode], strict=True):
            assert abs(rate - reference) <= 5.0 * math.hypot(error, errors["F"][0])


# Validate serialized VEGAS arrays and their diagnostic plots
def test_vgrid():
    ## Execute MC generator
    print(__name__ + ": Generating events ...")

    cmd = "./bin/gr -i gencard/test.json -p 'GP[CON]<F> -> pi+ pi-' -w true -l false -h 0 -n 0 -g VEGAS"
    execute(cmd)

    ## Read integration output
    with open("./vgrid/test.vgrid") as f:
        data = json.load(f)

    print(data.keys())
    print(data["VEGAS"].keys())

    vegas = data["VEGAS"]
    with open("./modeldata/TUNE0/NUMERICS.json", encoding="utf-8") as stream:
        configured_vegas = json5.load(stream)["NUMERICS_VEGAS"]

    assert data["SCHEMA_VERSION"] == 1
    assert data["PROPOSAL_ONLY"] is False
    for key in ("uniform_mix", "ess_rel_tolerance"):
        assert vegas[key] == pytest.approx(float(configured_vegas[key]))
    for key in (
        "min_support",
        "max_support_calls",
        "max_rounds",
        "convergence_window",
    ):
        assert vegas[key] == int(configured_vegas[key])
    assert vegas["best_grid_choice"] == configured_vegas["best_grid_choice"]
    assert (
        vegas["automatic_convergence"]
        is configured_vegas["automatic_convergence"]
    )
    assert "grid_rel_tolerance" not in vegas
    assert "min_bins" not in vegas
    assert "min_bin_width" not in vegas
    assert data["MAX_WEIGHT"]["max_w_value"] is None
    assert "mode" not in data["MAX_WEIGHT"]
    removed_keys = {
        'dcache', 'dxvec', 'region', 'rvec', 'xcache',
        'batch_mean', 'batch_m2', 'batch_calls', 'batch_finite', 'batch_committed',
        'estimate_mean', 'estimate_m2', 'estimate_weight', 'estimate_max_log_weight', 'estimate_batches',
        'continuation_base_mean', 'continuation_base_m2', 'continuation_base_weight',
        'continuation_base_max_log_weight', 'continuation_base_batches',
        'continuation_mean', 'continuation_m2', 'continuation_calls', 'continuation_finite', 'continuation_calibrated',
    }
    assert removed_keys.isdisjoint(vegas), removed_keys.intersection(vegas)

    statistics = data["STAT"]
    moments = ("integrand", "proposal", "importance", "cross_section")
    fields = ("count", "positive_count", "scaled_mean", "scaled_m2", "log_scale", "min_log_value", "finite")
    required = {f"{moment}_{field}" for moment in moments for field in fields}
    required.update((
        "integration_runtime", "integration_batch_mean", "integration_batch_m2",
        "chi2_estimate_mean", "chi2_estimate_scaled_m2", "chi2_estimate_scaled_weight",
    ))
    assert required <= statistics.keys(), required - statistics.keys()
    for moment in moments:
        expected_count = 0 if moment == "cross_section" else statistics["integration_samples"]
        assert statistics[f"{moment}_count"] == expected_count, moment
        assert statistics[f"{moment}_positive_count"] <= expected_count, moment
        assert statistics[f"{moment}_finite"] is True, moment
    assert statistics["proposal_positive_count"] == statistics["proposal_count"]
    assert statistics["cross_section_positive_count"] == 0
    removed_raw_statistics = {
        'Wsum', 'W2sum', 'sampling_Wsum', 'sampling_W2sum', 'sampling_count', 'sampling_scale',
        'sampling_scaled_sum', 'sampling_scaled_square_sum', 'sampling_scaled_absolute_sum', 'sampling_finite',
        'direct_count', 'direct_mean', 'direct_m2', 'direct_finite', 'sigma_err2', 'maxW', 'maxf',
        'chi2_weighted_sum', 'chi2_weighted_square_sum',
    }
    assert removed_raw_statistics.isdisjoint(statistics), removed_raw_statistics.intersection(statistics)

    diagnostic_dir = OUTPUT.plots / "vgrid"
    diagnostic_dir.mkdir(parents=True, exist_ok=True)

    keys = ["xmat", "f2mat", "fmat"]

    for key in keys:
        values = np.asarray(data["VEGAS"][key], dtype=float)
        assert values.size > 0
        assert np.isfinite(values).all()
        fig, ax = plt.subplots()
        plt.title(f"{key}")
        plt.plot(values)
        plot_path = diagnostic_dir / f"testbench_vgrid_{key}.pdf"
        fig.savefig(plot_path, bbox_inches="tight")
        plt.close(fig)
        assert plot_path.is_file()
        assert plot_path.stat().st_size > 0
