# Shared helpers for iceplot pytest workflows
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import codecs
import errno
import json
import math
import os
import pathlib
import pty
import select
import shlex
import signal
import subprocess
import sys
import time

import numpy as np

CDIR = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(CDIR / "python"))

from core.io import readers
from core.plot import plot as plots

from tests.technical.support.environment import project_environment_command
from tests.technical.support.output import PhysicsOutput

inputcard = str(CDIR / "icepack" / "SOFTCEP" / "processes" / "_common" / "gencard.json")
N_EVENTS = 3000
PROCESS = "GP[CON]<F> -> pi+ pi-"
PID = [[211, -211]]


# Stream a pseudo-terminal backed child process while keeping output for assertions
def stream_process_output(proc, master_fd, timeout):
    decoder = codecs.getincrementaldecoder("utf-8")("replace")
    chunks = []
    deadline = None if timeout is None else time.monotonic() + timeout

    while True:
        if deadline is not None and time.monotonic() > deadline:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            raise subprocess.TimeoutExpired(proc.args, timeout, output="".join(chunks))

        ready, _, _ = select.select([master_fd], [], [], 0.1)
        if master_fd in ready:
            try:
                data = os.read(master_fd, 4096)
            except OSError as exc:
                if exc.errno == errno.EIO:
                    break
                raise
            if not data:
                break
            text = decoder.decode(data)
            chunks.append(text)
            sys.stdout.write(text)
            sys.stdout.flush()
            continue

        if proc.poll() is not None:
            while True:
                ready, _, _ = select.select([master_fd], [], [], 0.0)
                if master_fd not in ready:
                    break
                try:
                    data = os.read(master_fd, 4096)
                except OSError as exc:
                    if exc.errno == errno.EIO:
                        break
                    raise
                if not data:
                    break
                text = decoder.decode(data)
                chunks.append(text)
                sys.stdout.write(text)
                sys.stdout.flush()
            break

    tail = decoder.decode(b"", final=True)
    if tail:
        chunks.append(tail)
        sys.stdout.write(tail)
        sys.stdout.flush()
    proc.wait()
    return "".join(chunks)


# Execute a command, stream output live and keep it in assertion failures
def run_command(cmd, timeout=1800):
    master_fd, slave_fd = pty.openpty()
    proc = None
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=CDIR,
            stdin=subprocess.DEVNULL,
            stdout=slave_fd,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            close_fds=True,
        )
        os.close(slave_fd)
        slave_fd = None
        output = stream_process_output(proc, master_fd, timeout)
    except subprocess.TimeoutExpired as exc:
        output = exc.output or ""
        raise AssertionError(f"Command timed out after {timeout} seconds\n{output}") from exc
    finally:
        if slave_fd is not None:
            os.close(slave_fd)
        os.close(master_fd)

    assert proc is not None
    assert proc.returncode == 0, output
    return output


# Run the graniitti generator through the project environment
def run_graniitti(args):
    quoted_args = " ".join(shlex.quote(arg) for arg in args)
    cmd = project_environment_command(f"./bin/gr {quoted_args}")
    return run_command(cmd)


# Run iceplot through the same project environment used by the generator
def run_iceplot(
    hepmc3_tags,
    plot_tag,
    labels,
    obs_module="default",
    cuts=None,
    pid=PID,
    unit="ub",
    analysis=None,
    mc_hepdata_samples=None,
    mc_scales=None,
    density=False,
    output=None,
):
    output = PhysicsOutput("iceplot") if output is None else output
    if not isinstance(output, PhysicsOutput):
        raise TypeError("run_iceplot output must be a PhysicsOutput")
    sources = [str(CDIR / "output" / f"{tag}.hepmc3") for tag in hepmc3_tags]
    pid_text = str(pid).replace(" ", "")
    cuts_arg = f"--cuts {' '.join(shlex.quote(cut) for cut in cuts)} " if cuts is not None else ""
    analysis_arg = f"--analysis {shlex.quote(analysis)} " if analysis is not None else ""
    mc_hepdata_sample_arg = (
        f"--mc-hepdata-sample {' '.join(shlex.quote(name) for name in mc_hepdata_samples)} "
        if mc_hepdata_samples is not None
        else ""
    )
    mc_scale_arg = (
        f"--mcscale {' '.join(shlex.quote(str(scale)) for scale in mc_scales)} "
        if mc_scales is not None
        else ""
    )
    density_arg = "--density " if density else ""
    command = (
        "python -m core.iceplot "
        f"--hepmc3 {' '.join(shlex.quote(source) for source in sources)} "
        f"--obs {shlex.quote(obs_module)} "
        f"{cuts_arg}"
        f"--pid {shlex.quote(pid_text)} "
        f"--mclabel {' '.join(shlex.quote(label) for label in labels)} "
        f"{analysis_arg}"
        f"{mc_hepdata_sample_arg}"
        f"{mc_scale_arg}"
        f"{density_arg}"
        f"--unit {shlex.quote(unit)} "
        f"--output {shlex.quote(plot_tag)} "
        f"--output-dir {shlex.quote(str(output.plots))} "
        f"--report {shlex.quote(str(output.report(plot_tag)))}"
    )
    cmd = project_environment_command(command)
    result = run_command(cmd)
    assert_iceplot_report(output.report(plot_tag))
    return result


# Load one iceplot comparison report
def load_iceplot_report(path):
    with pathlib.Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


# Require every reported sample and generated plot to be numerically valid
def assert_iceplot_report(path):
    report_path = pathlib.Path(path)
    assert report_path.is_file(), f"Missing iceplot report {report_path}"
    text_path = report_path.with_suffix(".txt")
    assert text_path.is_file(), f"Missing iceplot table {text_path}"
    report = load_iceplot_report(report_path)
    assert report.get("schema_version") == 1
    assert report.get("normalization") in {"cross_section", "unit_density"}
    assert report.get("ratio_uncertainty") in {"combined", "separate", "numerator", "none"}
    assert isinstance(report.get("stack"), bool)
    assert report.get("sets")
    plot_root = pathlib.Path(report["plot_directory"])
    populated_histograms = 0
    histogram_samples = 0
    for dataset in report["sets"]:
        assert {"title", "region", "mc_scale"}.issubset(dataset)
        assert dataset["title"] is None or isinstance(dataset["title"], str)
        assert dataset["region"] is None or isinstance(dataset["region"], str)
        assert math.isfinite(dataset["mc_scale"])
        assert dataset["observables"]
        dataset_root = plot_root / dataset["plotname"]
        for observable in dataset["observables"]:
            if observable["kind"] == "roc":
                plot = dataset_root / "roc" / f"rocplot__{observable['observable']}.pdf"
                assert plot.is_file(), f"Missing iceplot plot {plot}"
            else:
                for scale in ("linear", "log"):
                    plot = dataset_root / scale / f"hplot__{observable['observable']}.pdf"
                    assert plot.is_file(), f"Missing iceplot plot {plot}"
            for sample in observable["samples"]:
                if observable["kind"] == "roc":
                    assert sample["events"] > 0
                    assert math.isfinite(sample["weight_sum"])
                    continue
                histogram_samples += 1
                assert sample["integral"] >= 0.0
                assert sample["integral_error"] >= 0.0
                assert math.isfinite(sample["integral"])
                assert math.isfinite(sample["integral_error"])
                assert sample["empty"] == (sample["integral"] == 0.0)
                assert {
                    "mc_max_rel_uncertainty",
                    "mc_zero_bin_count",
                    "mc_valid_bin_count",
                }.issubset(sample)
                valid_bins = sample["mc_valid_bin_count"]
                zero_bins = sample["mc_zero_bin_count"]
                assert isinstance(valid_bins, int) and valid_bins >= 0
                assert isinstance(zero_bins, int) and 0 <= zero_bins <= valid_bins
                max_rel_uncertainty = sample["mc_max_rel_uncertainty"]
                if valid_bins == zero_bins:
                    assert max_rel_uncertainty is None
                else:
                    assert max_rel_uncertainty is not None
                    assert math.isfinite(max_rel_uncertainty)
                    assert max_rel_uncertainty >= 0.0
                source_fields = {
                    "mc_selected_events",
                    "mc_ess_fraction",
                    "mc_effective_events",
                }
                present_source_fields = source_fields.intersection(sample)
                assert not present_source_fields or present_source_fields == source_fields
                if present_source_fields:
                    selected_events = sample["mc_selected_events"]
                    ess_fraction = sample["mc_ess_fraction"]
                    effective_events = sample["mc_effective_events"]
                    assert isinstance(selected_events, int) and selected_events >= 0
                    assert math.isfinite(ess_fraction) and 0.0 <= ess_fraction <= 1.0
                    assert math.isfinite(effective_events) and effective_events >= 0.0
                    assert math.isclose(
                        effective_events,
                        ess_fraction * selected_events,
                        rel_tol=1.0e-12,
                        abs_tol=1.0e-12,
                    )
                populated_histograms += int(not sample["empty"])
                if sample.get("chi2") is not None:
                    assert sample["ndf"] > 0
                    assert math.isfinite(sample["chi2"])
                    assert math.isfinite(sample["chi2_ndf"])
                    if sample["empty"]:
                        assert sample["comparison_status"] == "empty_mc"
                        assert sample["shape_l1"] is None
                    else:
                        assert sample["comparison_status"] == "compared"
                        assert 0.0 <= sample["shape_l1"] <= 2.0
                elif "chi2" in sample:
                    assert sample["chi2_ndf"] is None
                    assert sample["chi2_status"] == "unavailable"
    assert histogram_samples == 0 or populated_histograms > 0, (
        "iceplot report contains no populated histograms"
    )
    return report


# Compute the cross-section scaled dPhi_pp histogram used to validate the overlay
def read_cross_section_histogram(tag):
    all_obs = readers.get_observables("default")
    obs = {"dPhi_pp": all_obs["dPhi_pp"].copy()}
    obs["dPhi_pp"]["bins"] = np.linspace(0.0, 180.0, 13)

    mcdata = readers.read_hepmc3(
        hepmc3file=str(CDIR / "output" / f"{tag}.hepmc3"),
        obs=[obs],
        pid=PID,
        cuts=["core.analysis.cuts.default"],
        chunk_range=[0, N_EVENTS - 1],
        verbose=0,
    )[0]

    return plots.histmc(
        mcdata=mcdata,
        obs=obs,
        density=False,
        scale=1e-6,
    )["dPhi_pp"]["hdata"]


# Compare cross-section histograms with relative integral and L1 distances
def assert_histograms_match(unweighted_hist, weighted_hist):
    np.testing.assert_allclose(unweighted_hist.bins, weighted_hist.bins, rtol=0.0, atol=1.0e-12)
    np.testing.assert_array_equal(unweighted_hist.valid, weighted_hist.valid)
    unweighted = unweighted_hist.counts_scaled
    weighted = weighted_hist.counts_scaled
    binwidth = unweighted_hist.binwidth

    assert np.isfinite(unweighted).all()
    assert np.isfinite(weighted).all()
    assert unweighted_hist.integral() > 0.0
    assert weighted_hist.integral() > 0.0

    errors = (unweighted_hist.integral_error(), weighted_hist.integral_error())
    for histogram, error in zip((unweighted_hist, weighted_hist), errors, strict=True):
        assert math.isfinite(error) and 0.0 < error < 0.05 * histogram.integral()
    assert abs(unweighted_hist.integral() - weighted_hist.integral()) <= 5.0 * math.hypot(*errors)

    mean_integral = 0.5 * (unweighted_hist.integral() + weighted_hist.integral())
    integral_distance = abs(unweighted_hist.integral() - weighted_hist.integral()) / mean_integral
    l1_distance = np.sum(np.abs(unweighted - weighted) * binwidth) / mean_integral

    assert integral_distance < 0.10
    assert l1_distance < 0.40
