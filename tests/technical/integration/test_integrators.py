# Monte Carlo integrator and event generation tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

# Run with: pytest -q -s -rP tests/technical/integration/test_integrators.py --run-integration

import json
import os
import pathlib
import shlex
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pyjson5 as json5
import pytest

from tests.technical.support.commands import execute
from tests.technical.support.environment import project_environment_command
from tests.technical.support.hepmc import read_events
from tests.technical.support.iceplot import run_iceplot
from tests.technical.support.output import PhysicsOutput

N = 400  # number of events
XS_UNIT = 1e6  # microbarns
GENCARD = "gencard/test.json"
OUTPUT = PhysicsOutput("generator_integrators")


# Read all event-level HepMC3 GenCrossSection values
def read_event_cross_sections(path):
    return read_event_cross_section_records(path)[:, 0]


# Read event-level HepMC3 cross sections and their sampling uncertainties
def read_event_cross_section_records(path):
    return np.asarray([(event.cross_section().xsec(), event.cross_section().xsec_err())
                       for event in read_events(path)])


# Compare generator completion and outputs across integration modes
def test_integrators():
    """
    Test MC integrators
    """

    ## Execute MC generator
    weightmode = ["true", "false"]
    integrator = ["VEGAS", "NEUROJAC"]

    for iw, w in enumerate(weightmode):
        for ig, g in enumerate(integrator):
            # Separate proposals and generated samples across all integration modes
            seed = 27183 + 7919 * (len(integrator) * iw + ig)
            seed_option = f" --set GENERIC.RNDSEED={seed}"
            reload_seed_option = f" --set GENERIC.RNDSEED={seed + 1000003}"
            print(f"Generating events with weightmode <{w}> and integrator <{g}> ...")
            neurojac_options = ""
            if g == "NEUROJAC":
                neurojac_options = (
                    " --set NUMERICS.json:NUMERICS_NEUROJAC.vegas_ncall=5000"
                    " --set NUMERICS.json:NUMERICS_NEUROJAC.vegas_rounds=2"
                    " --set NUMERICS.json:NUMERICS_NEUROJAC.vegas_sampler_rounds=1"
                )

            # Generate and reload each serialized importance proposal
            cmd = f"./bin/gr -i {GENCARD} -p 'GP[CON]<F> -> pi+ pi-' -l false -h 0 -c 4 -w {w} -n 0 -g {g} -o {g}_w_{w}{neurojac_options}"
            execute(cmd + seed_option)
            cmd = f"./bin/gr -i {GENCARD} -p 'GP[CON]<F> -> pi+ pi-' -l false -h 0 -c 4 -w {w} -n {N} -g {g} -o {g}_w_{w} -d vgrid/{g}_w_{w}.vgrid"
            execute(cmd + reload_seed_option)

    # Make plots
    run_iceplot(
        hepmc3_tags=[
            "VEGAS_w_false",
            "VEGAS_w_true",
            "NEUROJAC_w_false",
            "NEUROJAC_w_true",
        ],
        plot_tag="testbench_integrators",
        labels=[
            "VEGAS unweighted",
            "VEGAS weighted",
            "NEUROJAC unweighted",
            "NEUROJAC weighted",
        ],
        pid=[[211, -211]] * 4,
        unit="ub",
        output=OUTPUT,
    )

    generated_tags = (
        "VEGAS_w_false",
        "VEGAS_w_true",
        "NEUROJAC_w_false",
        "NEUROJAC_w_true",
    )
    terminal_estimates = {}
    card = json5.loads(pathlib.Path(GENCARD).read_text(encoding="utf-8"))
    for tag in generated_tags:
        event_file = pathlib.Path("output") / f"{tag}.hepmc3"
        assert event_file.is_file()
        assert event_file.stat().st_size > 0
        estimates = read_event_cross_section_records(event_file)
        assert estimates.shape == (N, 2)
        assert np.isfinite(estimates).all()
        assert np.all(estimates > 0.0)
        assert estimates[-1, 0] / estimates[0, 0] == pytest.approx(1.0, rel=0.10)
        # A sample variance can increase when a rare large weight is observed
        precision = card["INTEGRATOR"]["precision"]
        assert estimates[-1, 1] / estimates[-1, 0] <= precision
        terminal_estimates[tag] = estimates[-1]

    for mode in integrator:
        weighted = terminal_estimates[f"{mode}_w_true"]
        unweighted = terminal_estimates[f"{mode}_w_false"]
        combined_error = np.hypot(weighted[1], unweighted[1])
        assert abs(weighted[0] - unweighted[0]) <= 5.0 * combined_error

    reference = terminal_estimates["VEGAS_w_true"]
    for tag, estimate in terminal_estimates.items():
        combined_error = np.hypot(reference[1], estimate[1])
        assert abs(reference[0] - estimate[0]) <= 5.0 * combined_error, tag

    for mode in ("VEGAS", "NEUROJAC"):
        for weight in weightmode:
            grid_file = pathlib.Path("vgrid") / f"{mode}_w_{weight}.vgrid"
            with grid_file.open(encoding="utf-8") as stream:
                grid = json.load(stream)
            assert grid["STAT"]["generated"] == N
            assert grid["STAT"]["trials"] >= N
            assert grid["STAT"]["integrand_positive_count"] > 0
            if mode == "NEUROJAC":
                assert grid["STAT"]["importance_positive_count"] > 0
            assert (
                grid["STAT"]["cross_section_count"]
                > grid["STAT"]["importance_count"]
            )

    plot_dir = OUTPUT.plot_dir("testbench_integrators")
    assert list(plot_dir.glob("*.pdf"))


# Check proposal-only mode saves VEGAS and NEUROJAC without production integration
@pytest.mark.parametrize(
    ("integrator", "overrides"),
    (
        ("VEGAS", ["NUMERICS.json:NUMERICS_VEGAS.automatic_convergence=false"]),
        (
            "NEUROJAC",
            [
                "NUMERICS.json:NUMERICS_NEUROJAC.buffer_size=256",
                "NUMERICS.json:NUMERICS_NEUROJAC.batch_size=128",
                "NUMERICS.json:NUMERICS_NEUROJAC.val_size=128",
                "NUMERICS.json:NUMERICS_NEUROJAC.rounds=1",
                "NUMERICS.json:NUMERICS_NEUROJAC.epochs=1",
                "NUMERICS.json:NUMERICS_NEUROJAC.vegas_ncall=5000",
                "NUMERICS.json:NUMERICS_NEUROJAC.vegas_rounds=1",
                "NUMERICS.json:NUMERICS_NEUROJAC.vegas_sampler_rounds=1",
            ],
        ),
    ),
)
def test_proposal_only_mode_saves_grid_without_xs(integrator, overrides):
    output_tag = f"test_{integrator.lower()}_proposal_only_{os.getpid()}"
    overrides = [
        *overrides,
        "INTEGRATOR.min_samples=2048",
        "INTEGRATOR.max_samples=8192",
        "INTEGRATOR.precision=0.9",
    ]
    arguments = [
        "./bin/gr",
        "-i",
        GENCARD,
        "-n",
        "-1",
        "-g",
        integrator,
        "-p",
        "GP[CON]<F> -> pi+ pi-",
        "-l",
        "false",
        "-w",
        "true",
        "-c",
        "2",
        "-h",
        "0",
        "-o",
        output_tag,
    ]
    for override in overrides:
        arguments.extend(("--set", override))
    result = subprocess.run(
        project_environment_command(" ".join(shlex.quote(argument) for argument in arguments)),
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=180,
    )

    assert result.returncode == 0, result.stdout
    with open(f"vgrid/{output_tag}.vgrid", encoding="utf-8") as stream:
        grid = json.load(stream)
    assert grid["STAT"]["integration_samples"] == 0
    assert grid["STAT"]["integration_runtime"] == pytest.approx(0.0)
    assert grid["STAT"]["sigma"] == pytest.approx(0.0)
    assert grid["MAX_WEIGHT"]["initialized"] is False
    assert grid["PROPOSAL_ONLY"] is True
    # Weighted production uses the proposal directly without hidden integration
    generation_arguments = list(arguments)
    generation_arguments[generation_arguments.index("-n") + 1] = "2"
    generation_arguments[generation_arguments.index("-l") + 1] = "true"
    generation_arguments[generation_arguments.index("-o") + 1] = f"{output_tag}_events"
    generation_arguments.extend(("-d", f"vgrid/{output_tag}.vgrid"))
    generation = subprocess.run(
        project_environment_command(
            " ".join(shlex.quote(argument) for argument in generation_arguments)
        ),
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=180,
    )
    assert generation.returncode == 0, generation.stdout
    event_path = pathlib.Path("output") / f"{output_tag}_events.hepmc3"
    event_cross_sections = read_event_cross_section_records(event_path)
    assert event_cross_sections.shape == (2, 2)
    assert np.isfinite(event_cross_sections[:, 0]).all()
    assert np.all(event_cross_sections[:, 0] > 0.0)
    assert np.isfinite(event_cross_sections[-1]).all()
    assert event_cross_sections[-1, 1] >= 0.0

    # Generation adds cross section samples without changing proposal-only state
    with open(f"vgrid/{output_tag}_events.vgrid", encoding="utf-8") as stream:
        generation_grid = json.load(stream)
    assert generation_grid["STAT"]["integration_samples"] == 0
    assert generation_grid["STAT"]["cross_section_count"] >= 2
    assert generation_grid["STAT"]["sigma"] > 0.0
    assert np.isfinite(generation_grid["STAT"]["sigma_err"])
    assert generation_grid["PROPOSAL_ONLY"] is True

    # An explicit zero-event reload completes and persists the missing integral
    integration_arguments = list(arguments)
    integration_arguments[integration_arguments.index("-n") + 1] = "0"
    integration_arguments.extend(("-d", f"vgrid/{output_tag}.vgrid"))
    integration = subprocess.run(
        project_environment_command(
            " ".join(shlex.quote(argument) for argument in integration_arguments)
        ),
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=180,
    )
    assert integration.returncode == 0, integration.stdout
    with open(f"vgrid/{output_tag}.vgrid", encoding="utf-8") as stream:
        integrated_grid = json.load(stream)
    assert integrated_grid["STAT"]["integration_samples"] >= 2048
    assert integrated_grid["STAT"]["integration_runtime"] > 0.0
    assert integrated_grid["STAT"]["sigma"] > 0.0
    assert integrated_grid["PROPOSAL_ONLY"] is False


# Check full integration retains recoverable cross section and maximum weight data
def test_zero_event_mode_full_integration_state():
    output_tag = f"test_vegas_full_integration_{os.getpid()}"
    overrides = [
        "INTEGRATOR.VEGAS.rounds=1",
        "INTEGRATOR.VEGAS.ncall=5000",
        "INTEGRATOR.min_samples=2048",
        "INTEGRATOR.max_samples=8192",
        "INTEGRATOR.precision=0.9",
        "NUMERICS.json:NUMERICS_VEGAS.automatic_convergence=false",
        "NUMERICS.json:NUMERICS_MC.MAX_WEIGHT.max_overflow_prob=0.5",
        "NUMERICS.json:NUMERICS_MC.MAX_WEIGHT.confidence_level=0.75",
    ]
    arguments = [
        "./bin/gr",
        "-i",
        GENCARD,
        "-n",
        "0",
        "-g",
        "VEGAS",
        "-p",
        "GP[CON]<F> -> pi+ pi-",
        "-l",
        "false",
        "-w",
        "false",
        "-c",
        "2",
        "-h",
        "0",
        "-o",
        output_tag,
    ]
    for override in overrides:
        arguments.extend(("--set", override))
    command = project_environment_command(
        " ".join(shlex.quote(argument) for argument in arguments)
    )
    result = subprocess.run(
        command,
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=180,
    )

    assert result.returncode == 0, result.stdout
    with open(f"vgrid/{output_tag}.vgrid", encoding="utf-8") as stream:
        grid = json.load(stream)
    assert grid["STAT"]["integration_samples"] >= 2048
    assert grid["STAT"]["integration_runtime"] > 0.0
    assert grid["STAT"]["sigma"] > 0.0
    assert grid["STAT"]["sigma_err"] >= 0.0
    assert grid["MAX_WEIGHT"]["initialized"] is True
    assert grid["MAX_WEIGHT"]["envelope"] > 0.0
    assert grid["PROPOSAL_ONLY"] is False
    assert grid["STAT"]["cross_section_count"] == 0

    reload_arguments = arguments + ["-d", f"vgrid/{output_tag}.vgrid"]
    reload_result = subprocess.run(
        project_environment_command(
            " ".join(shlex.quote(argument) for argument in reload_arguments)
        ),
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=180,
    )
    assert reload_result.returncode == 0, reload_result.stdout

    mismatch_arguments = list(reload_arguments)
    mismatch_arguments[mismatch_arguments.index("-l") + 1] = "true"
    mismatch = subprocess.run(
        project_environment_command(
            " ".join(shlex.quote(argument) for argument in mismatch_arguments)
        ),
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=180,
    )
    assert mismatch.returncode != 0


# Check invalid proposal-only combinations fail before numerical sampling
@pytest.mark.parametrize(
    "arguments",
    (
        ("-n", "-2"),
        ("-n", "-1", "-g", "FLAT"),
        ("-n", "-1", "-d", "vgrid/missing.vgrid"),
    ),
)
def test_proposal_only_mode_invalid_cases(arguments):
    command = [
        "./bin/gr",
        "-i",
        GENCARD,
        "-p",
        "GP[CON]<F> -> pi+ pi-",
        "-l",
        "false",
        "-h",
        "0",
        *arguments,
    ]
    result = subprocess.run(
        project_environment_command(
            " ".join(shlex.quote(argument) for argument in command)
        ),
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=120,
    )

    assert result.returncode != 0


# Reject an unknown integration algorithm
def test_invalid_integrator_is_rejected():
    arguments = [
        "./bin/gr",
        "-i",
        GENCARD,
        "-n",
        "0",
        "-g",
        "NEUROJA",
        "-h",
        "0",
    ]
    result = subprocess.run(
        project_environment_command(" ".join(shlex.quote(argument) for argument in arguments)),
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=120,
    )

    assert result.returncode != 0


# Require an impossible phase space to fail before accepting VEGAS adaptation
def test_vegas_zero_support_is_fatal():
    output_tag = f"test_vegas_zero_support_{os.getpid()}"
    arguments = [
        "./bin/gr",
        "-i",
        GENCARD,
        "-n",
        "0",
        "-g",
        "VEGAS",
        "-e",
        "1",
        "-p",
        "GP[CON]<F> -> pi+ pi-",
        "-o",
        output_tag,
        "-h",
        "0",
    ]
    command = project_environment_command(" ".join(shlex.quote(argument) for argument in arguments))
    result = subprocess.run(
        command,
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=300,
    )

    assert result.returncode != 0, result.stdout
    assert not (pathlib.Path("vgrid") / f"{output_tag}.vgrid").exists()


# Require direct and learned proposals to fail explicitly on zero support
@pytest.mark.parametrize(
    ("integrator", "overrides"),
    [
        (
            "FLAT",
            ["INTEGRATOR.min_samples=1000"],
        ),
        (
            "NEUROJAC",
            [
                "NUMERICS.json:NUMERICS_NEUROJAC.buffer_size=1",
                "NUMERICS.json:NUMERICS_NEUROJAC.vegas_ncall=50",
                "NUMERICS.json:NUMERICS_NEUROJAC.vegas_rounds=1",
                "NUMERICS.json:NUMERICS_NEUROJAC.vegas_sampler_rounds=1",
            ],
        ),
    ],
)
def test_other_integrators_zero_support_fatal(integrator, overrides):
    output_tag = f"test_{integrator.lower()}_zero_support_{os.getpid()}"
    arguments = [
        "./bin/gr",
        "-i",
        GENCARD,
        "-n",
        "0",
        "-e",
        "1",
        "-p",
        "GP[CON]<F> -> pi+ pi-",
        "-g",
        integrator,
        "-o",
        output_tag,
        "-h",
        "0",
    ]
    for override in overrides:
        arguments.extend(("--set", override))
    command = project_environment_command(" ".join(shlex.quote(argument) for argument in arguments))
    result = subprocess.run(
        command,
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=120,
    )
    assert result.returncode != 0, result.stdout
    assert not (pathlib.Path("vgrid") / f"{output_tag}.vgrid").exists()


# Require forward-veto trials to remain event-local and improve the VEGAS estimate
def test_vegas_forward_veto_generation_continuation():
    output_tag = f"test_vegas_forward_veto_{os.getpid()}"
    event_path = pathlib.Path("output") / f"{output_tag}.hepmc3"
    grid_path = pathlib.Path("vgrid") / f"{output_tag}.vgrid"
    # Preserve the same cuts and integration settings when reloading the grid
    settings = [
        '--set', 'SCATTERING.BEAMFRAG="cylinder"',
        "--set", "VETOCUTS.active=true",
        "--set", "INTEGRATOR.VEGAS.ncall=30000",
        "--set", "INTEGRATOR.min_samples=30000",
        "--set", "INTEGRATOR.precision=0.1",
    ]
    arguments = [
        "./bin/gr",
        "-i",
        "icepack/GAMMA/ATLAS_1377585/mumu/gencard.json",
        "-n",
        "400",
        "-g",
        "VEGAS",
        "-w",
        "true",
        "-s",
        "1",
        "-c",
        "4",
        "-h",
        "0",
        "-o",
        output_tag,
        *settings,
    ]
    command = project_environment_command(" ".join(shlex.quote(argument) for argument in arguments))
    result = subprocess.run(
        command,
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        # Full VEGAS adaptation includes forward fragmentation for the veto
        timeout=3600,
    )
    assert result.returncode == 0, result.stdout

    values = read_event_cross_sections(event_path)
    assert values.size == 400
    assert np.isfinite(values).all()
    assert np.all(values > 0.0)
    assert values[-1] / values[0] == pytest.approx(1.0, rel=0.10)

    with grid_path.open(encoding="utf-8") as stream:
        grid = json.load(stream)
    assert grid["STAT"]["generated"] == 400
    assert grid["STAT"]["trials"] >= 400
    assert (
        grid["STAT"]["cross_section_count"]
        > grid["STAT"]["importance_count"]
    )
    veto_efficiency = grid["STAT"]["vetocuts_ok"] / grid["STAT"]["fidcuts_ok"]
    assert 0.90 < veto_efficiency < 1.0

    reload_tag = f"{output_tag}_reload"
    reload_arguments = [
        "./bin/gr",
        "-i",
        "icepack/GAMMA/ATLAS_1377585/mumu/gencard.json",
        "-n",
        "20",
        "-g",
        "VEGAS",
        "-w",
        "true",
        "-s",
        "1",
        "-c",
        "4",
        "-h",
        "0",
        "-d",
        str(grid_path),
        "-o",
        reload_tag,
        *settings,
    ]
    reload_command = project_environment_command(
        " ".join(shlex.quote(argument) for argument in reload_arguments)
    )
    reload_result = subprocess.run(
        reload_command,
        cwd=pathlib.Path(__file__).resolve().parents[3],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=120,
    )
    assert reload_result.returncode == 0, reload_result.stdout
    reload_values = read_event_cross_sections(pathlib.Path("output") / f"{reload_tag}.hepmc3")
    assert reload_values.size == 20
    assert reload_values[-1] / values[-1] == pytest.approx(1.0, rel=0.10)


# Compare independent VEGAS differential cross sections using stored squared weights
def test_hfast():
    spectra = []
    fig, ax = plt.subplots()
    for seed in (27183, 35102):
        tag = f"test_hfast_{os.getpid()}_{seed}"
        execute(f"./bin/gr -i {GENCARD} -p 'GP[CON]<F> -> pi+ pi-' -w true -l false -h 1 -n 0 -c 1 -g VEGAS -r {seed} -o {tag}")
        data = json.loads(pathlib.Path(f"output/{tag}.hfast").read_text())["h1"]["dPhi_pp"]
        edges = np.asarray(data["binedges"])
        widths = edges[:, 1] - edges[:, 0]
        weights = np.asarray(data["weights"])
        squares = np.asarray(data["weights2"])
        count = data["fills"]
        assert count > 1 and np.all(widths > 0.0)
        assert np.isfinite(weights).all() and np.isfinite(squares).all()
        assert np.all(weights >= 0.0) and np.any(weights > 0.0)
        assert np.all(squares >= 0.0)
        integral = weights / count * XS_UNIT
        stat = json.loads(pathlib.Path(f"vgrid/{tag}.vgrid").read_text())["STAT"]
        # Histogram accumulation starts after the common bin ranges are initialized
        assert count <= stat["integration_samples"]
        # Poisson sum-of-squares errors conservatively include the finite trial covariance
        variance = squares / count**2 * XS_UNIT**2
        assert np.all(variance[weights > 0.0] > 0.0)
        error = np.hypot(np.sqrt(variance.sum()), stat["sigma_err"] * XS_UNIT)
        assert abs(integral.sum() - stat["sigma"] * XS_UNIT) <= 5.0 * error
        spectra.append((edges, integral, variance))
        ax.errorbar(edges.mean(axis=1), integral / widths, yerr=np.sqrt(variance) / widths,
                    label=f"VEGAS seed {seed}", fmt=".")
    np.testing.assert_allclose(spectra[0][0], spectra[1][0])
    delta = spectra[0][1] - spectra[1][1]
    variance = spectra[0][2] + spectra[1][2]
    assert np.all(np.abs(delta) <= 5.0 * np.sqrt(variance))
    assert abs(delta.sum()) <= 5.0 * np.sqrt(variance.sum())
    ax.set(xlabel=r"$|\Delta\phi_{pp}|$ [deg]", ylabel=r"$d\sigma/d|\Delta\phi_{pp}|$ [$\mu$b/deg]")
    ax.legend()
    hfast_dir = OUTPUT.plots / "hfast"
    hfast_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(hfast_dir / "testbench_hfast.pdf", bbox_inches="tight")
    plt.close(fig)
