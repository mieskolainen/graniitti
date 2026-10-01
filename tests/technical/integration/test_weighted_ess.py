# Tests for weighted effective sample size output
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os
import pathlib
import shlex
import subprocess

import numpy as np

from tests.technical.support.environment import project_environment_command
from tests.technical.support.hepmc import read_events

CDIR = pathlib.Path(__file__).resolve().parents[3]
INPUTCARD = str(CDIR / "gencard" / "test.json")
PROCESS = "GP[CON]<F> -> pi+ pi-"
N_EVENTS = 32


# Run graniitti through the active or resolved project conda environment
def run_graniitti_conda(args):
    quoted_args = " ".join(shlex.quote(arg) for arg in args)
    cmd = project_environment_command(f"./bin/gr {quoted_args}")
    proc = subprocess.run(
        cmd,
        cwd=CDIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=900,
    )
    assert proc.returncode == 0, proc.stdout
    return proc.stdout


# Read generated weights and cumulative ESS through the official HepMC3 binding
def read_weighted_ess_payload(path):
    weights = []
    ess_values = []
    ess_fraction_values = []

    for event in read_events(path):
        assert len(event.weights()) == 1
        weights.append(event.weights()[0])
        ess_values.append(float(event.attribute_as_string("weighted_ess")))
        ess_fraction_values.append(float(event.attribute_as_string("weighted_ess_fraction")))

    return weights, ess_values, ess_fraction_values


# Validate every cumulative ESS against the generated importance weights
def test_run_reports_effective_sample_size():
    suffix = str(os.getpid())
    grid_tag = f"test_weighted_ess_grid_{suffix}"
    weighted_tag = f"test_weighted_ess_{suffix}"

    common_args = [
        "-i",
        INPUTCARD,
        "-p",
        PROCESS,
        "-l",
        "false",
        "-h",
        "0",
        "-f",
        "hepmc3",
        "-c",
        "1",
        "-g",
        "VEGAS",
        "--set", "INTEGRATOR.min_samples=1000",
        "--set", "INTEGRATOR.max_samples=2000",
        "--set", "INTEGRATOR.precision=1.0",
        "--set", "INTEGRATOR.VEGAS.ncall=1000",
        "--set", "INTEGRATOR.VEGAS.rounds=2",
        "--set", "GENCUTS.<F>.M=[0.5,2.0]",
    ]

    run_graniitti_conda(common_args + ["-w", "true", "-n", "0", "-o", grid_tag])
    run_graniitti_conda(
        common_args
        + ["-w", "true", "-n", str(N_EVENTS), "-d", f"vgrid/{grid_tag}.vgrid", "-o", weighted_tag]
    )

    weights, ess_values, ess_fraction_values = read_weighted_ess_payload(
        CDIR / "output" / f"{weighted_tag}.hepmc3"
    )
    assert len(weights) == N_EVENTS
    assert len(ess_values) == N_EVENTS
    assert len(ess_fraction_values) == N_EVENTS

    weights = np.asarray(weights)
    assert np.isfinite(weights).all() and (weights > 0.0).all()
    weights /= np.max(weights)
    assert not np.allclose(weights, weights[0], rtol=1e-12, atol=0.0)
    expected_ess = np.cumsum(weights)**2 / np.cumsum(weights**2)
    expected_fraction = expected_ess / np.arange(1, N_EVENTS + 1)
    np.testing.assert_allclose(ess_values, expected_ess, rtol=1e-12)
    np.testing.assert_allclose(ess_fraction_values, expected_fraction, rtol=1e-12)
