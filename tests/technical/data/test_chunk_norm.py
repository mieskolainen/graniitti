# Test iceplot chunking and cross-section normalization
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import gzip
import sys
import types
from functools import partial

import numpy as np
import pytest
from core.io import readers
from core.io.cache import cache

from tests.technical.support.hepmc import write_muon_events


# Share the official HepMC3 event writer across the chunk-normalization cases
@pytest.fixture
def hepmc_file(tmp_path):
    return partial(write_muon_events, tmp_path / 'events.hepmc3')


# Write a cut module consumed by the normal module loader
@pytest.fixture
def cut_file(tmp_path):
    # Serialize a controlled event selection and its parameters
    def write(name, parameters, body='return True'):
        path = tmp_path / f'{name}.py'
        path.write_text(
            '# Test event selection\n#\n# (c) 2026 Mikael Mieskolainen\n'
            '# Licensed under the MIT License <http://opensource.org/licenses/MIT>.\n\n'
            f'cut_param = {parameters!r}\n\n# Apply the event selection\n'
            'def cut_func(event):\n    ' + body.replace('\n', '\n    ') + '\n'
        )
        return str(path)
    return write


# Observe real C++ reader operations without substituting its event or reader classes
@pytest.fixture
def reader_calls(monkeypatch):
    calls = {"skip": [], "read": [], "close": []}
    reader = readers.hepmc3.ReaderAscii
    skip, read_event, close = reader.skip, reader.read_event, reader.close

    # Record the exact skip argument passed to the official reader
    def record_skip(self, count):
        calls["skip"].append(count)
        return skip(self, count)

    # Record successfully decoded event numbers
    def record_read(self, event):
        result = read_event(self, event)
        if not self.failed():
            calls["read"].append(event.event_number())
        return result

    # Record explicit closure on ordinary and exceptional exits
    def record_close(self):
        calls["close"].append(True)
        return close(self)

    monkeypatch.setattr(reader, "skip", record_skip)
    monkeypatch.setattr(reader, "read_event", record_read)
    monkeypatch.setattr(reader, "close", record_close)
    return calls


# Load the extensionless iceplot entry point as an importable test module
def load_iceplot_module():
    from importlib import import_module

    return import_module('core.iceplot')


# Exercise reconstruction failures on real serialized events
def run_reconstruction_failure_case(capsys, hepmc_file, verbose):
    # Compute a defined observable for every event
    def stable_observable(event):
        return event.evt.event_number()

    # Raise on a controlled subset with an undefined reconstructed observable
    def flaky_observable(event):
        if event.evt.event_number() < 2:
            raise ValueError(f"missing cascade for event {event.evt.event_number()}")
        return event.evt.event_number() + 10

    output = readers.read_hepmc3(
        hepmc3file=hepmc_file(),
        obs=[
            {
                "stable": {"func": stable_observable},
                "flaky": {"func": flaky_observable},
            }
        ],
        pid=[[211, -211]],
        cuts=[None],
        chunk_range=[0, 2],
        xsmode="header",
        verbose=verbose,
    )

    return output, capsys.readouterr().out.splitlines()


# Check different post-cuts share only explicitly cached event projections
def test_hepmc_cache_before_postcuts(cut_file, hepmc_file):
    calls = []

    # Record the physical events for which the shared projection is evaluated
    @cache
    def project(event):
        """Return the source event number."""
        calls.append(event.evt.event_number())
        return event.evt.event_number()

    selections = [
        cut_file('even', {'scale': 1.0}, 'return event.evt.event_number() % 2 == 0'),
        cut_file('above', {'scale': 1.0}, 'return event.evt.event_number() >= 1'),
    ]

    output = readers.read_hepmc3(
        hepmc3file=hepmc_file(),
        obs=[{"event": {"func": project}}, {"event": {"func": project}}],
        pid=[[211, -211], [211, -211]],
        cuts=selections,
        chunk_range=[0, 2],
        xsmode="header",
    )

    assert calls == [0, 1, 2]
    assert output[0]["data"]["event"].tolist() == [0, 2]
    assert output[1]["data"]["event"].tolist() == [1, 2]
    assert output[0]["event_ids"].tolist() == [0, 2]
    assert output[1]["event_ids"].tolist() == [1, 2]


# Check particle identities and cut parameters isolate event projector values
def test_hepmc_isolates_different_event_views(cut_file, hepmc_file):
    calls = []

    # Compute a value sensitive to every cache-group input
    @cache
    def project(event):
        """Return one PID and cut-parameter dependent event value."""
        context = (event.pid[0], event.cut_param["offset"])
        calls.append((event.evt.event_number(), *context))
        return event.evt.event_number() + context[0] + context[1]

    selections = [cut_file(name, {'offset': offset})
                  for name, offset in [('first', 1), ('second', 2), ('third', 1)]]

    output = readers.read_hepmc3(
        hepmc3file=hepmc_file(),
        obs=[{"value": {"func": project}}] * 3,
        pid=[[13, -13], [13, -13], [211, -211]],
        cuts=selections,
        chunk_range=[0, 2],
        xsmode="header",
    )

    assert len(calls) == 9
    assert output[0]["data"]["value"].tolist() == [14, 15, 16]
    assert output[1]["data"]["value"].tolist() == [15, 16, 17]
    assert output[2]["data"]["value"].tolist() == [212, 213, 214]


# Check shared selection parameters cannot mutate during post-cut evaluation
def test_hepmc_rejects_shared_param_mutation(cut_file, hepmc_file, reader_calls):
    selections = [
        cut_file('mutating', {'offset': 1}, "event.cut_param['offset'] = 2\nreturn True"),
        cut_file('stable', {'offset': 1}),
    ]

    with pytest.raises(RuntimeError, match="must not change"):
        readers.read_hepmc3(
            hepmc3file=hepmc_file(),
            obs=[{"event": {"func": lambda event: 1.0}}] * 2,
            pid=[[211, -211], [211, -211]],
            cuts=selections,
            chunk_range=[0, 2],
            xsmode="header",
        )
    assert reader_calls["close"] == [True]


# Check unsupported mutable parameter types conservatively disable cache sharing
def test_selection_cache_isolate_unsupported_params():
    selections = [
        types.SimpleNamespace(cut_param={"window": np.asarray([0.0, 1.0])}),
        types.SimpleNamespace(cut_param={"window": np.asarray([0.0, 1.0])}),
    ]

    cache_ids, cache_count, _ = readers._selection_cache_ids(
        pid=[[211, -211], [211, -211]],
        selections=selections,
    )

    assert cache_ids == [0, 1]
    assert cache_count == 2


# Check exact nested parameters preserve mapping order and signed zeros
def test_selection_cache_exact_structural_params():
    selections = [
        types.SimpleNamespace(cut_param={"active": True, "window": [0.0, 1.0]}),
        types.SimpleNamespace(cut_param={"window": [0.0, 1.0], "active": True}),
        types.SimpleNamespace(cut_param={"zero": -0.0}),
        types.SimpleNamespace(cut_param={"zero": 0.0}),
    ]

    cache_ids, cache_count, _ = readers._selection_cache_ids(
        pid=[[211, -211]] * 4,
        selections=selections,
    )

    assert cache_ids == [0, 1, 2, 3]
    assert cache_count == 4


# Check custom container state cannot enter a shared selection cache
def test_cache_container_subclasses():
    # Attach observable state outside the base dictionary contents
    class Parameters(dict):
        def __init__(self, label):
            super().__init__({"value": 1})
            self.label = label

    selections = [
        types.SimpleNamespace(cut_param=Parameters("first")),
        types.SimpleNamespace(cut_param=Parameters("second")),
    ]

    cache_ids, cache_count, _ = readers._selection_cache_ids(
        pid=[[211, -211], [211, -211]],
        selections=selections,
    )

    assert cache_ids == [0, 1]
    assert cache_count == 2


# Check non-finite parameters conservatively receive isolated caches
def test_selection_cache_isolate_nonfinite_params():
    selections = [
        types.SimpleNamespace(cut_param={"limit": np.nan}),
        types.SimpleNamespace(cut_param={"limit": np.nan}),
    ]

    cache_ids, cache_count, _ = readers._selection_cache_ids(
        pid=[[211, -211], [211, -211]],
        selections=selections,
    )

    assert cache_ids == [0, 1]
    assert cache_count == 2


# Assert inclusive chunk ranges cover every event exactly once
def assert_complete_partitions(partitions, totalsize):
    assert partitions[0][0] == 0
    assert partitions[-1][1] == totalsize - 1
    assert sum(end - start + 1 for start, end in partitions) == totalsize
    assert all(
        partitions[index][1] + 1 == partitions[index + 1][0] for index in range(len(partitions) - 1)
    )


# Check automatic chunks use the CPU limit without creating microsplits
def test_iceplot_auto_chunks_balanced_complete():
    module = load_iceplot_module()
    partitions, workers = module.build_event_partitions(
        nevents=9001,
        chunksize="auto",
        cores=64,
    )
    sizes = [end - start + 1 for start, end in partitions]

    assert workers == 9
    assert len(partitions) == 9
    assert min(sizes) == 1000
    assert max(sizes) == 1001
    assert_complete_partitions(partitions, totalsize=9001)


# Check the core count limits automatic parallel partitions
def test_iceplot_auto_chunks_respect_core_count():
    module = load_iceplot_module()
    partitions, workers = module.build_event_partitions(
        nevents=10500,
        chunksize="auto",
        cores=4,
    )

    assert workers == 4
    assert [end - start + 1 for start, end in partitions] == [2625] * 4
    assert_complete_partitions(partitions, totalsize=10500)


# Check a sample smaller than the minimum remains one complete chunk
def test_iceplot_auto_chunks_small_sample_complete():
    module = load_iceplot_module()
    partitions, workers = module.build_event_partitions(
        nevents=999,
        chunksize="auto",
        cores=16,
    )

    assert partitions == [(0, 998)]
    assert workers == 1
    assert_complete_partitions(partitions, totalsize=999)


# Check fixed chunking absorbs a sub-minimum trailing remainder
def test_iceplot_fixed_chunks_avoid_small_remainder():
    module = load_iceplot_module()
    partitions, workers = module.build_event_partitions(
        nevents=2500,
        chunksize=1000,
        cores=8,
    )

    assert partitions == [(0, 999), (1000, 2499)]
    assert workers == 2
    assert_complete_partitions(partitions, totalsize=2500)


# Check runtime coverage validation rejects missing or duplicated events
def test_noncontiguous_partitions():
    module = load_iceplot_module()

    with pytest.raises(RuntimeError, match="non-contiguous"):
        module.validate_event_partitions([(0, 999), (1001, 2000)], nevents=2001)
    with pytest.raises(RuntimeError, match="non-contiguous"):
        module.validate_event_partitions([(0, 999), (999, 1999)], nevents=2000)


# Check automatic chunking is the CLI default and explicit values are bounded
def test_iceplot_chunksize_cli_default_validation(monkeypatch, hepmc_file):
    module = load_iceplot_module()
    monkeypatch.setattr(sys, "argv", ["iceplot", "--hepmc3", hepmc_file()])

    args = module.parse_args()

    assert args.chunksize == "auto"
    assert module.parse_chunksize("1000") == 1000
    with pytest.raises(module.argparse.ArgumentTypeError, match="at least 1000"):
        module.parse_chunksize("999")


def test_trial_counter_difference():
    chunk_trials = readers._chunk_attempted_events(
        chunk_range=[1000, 1999],
        attempted_before_chunk=1353,
        attempted_in_chunk_last=2696,
    )

    assert chunk_trials == pytest.approx(1343)


def test_attempted_events_handles_full_range_zero():
    chunk_trials = readers._chunk_attempted_events(
        chunk_range=[0, 3999],
        attempted_before_chunk=0,
        attempted_in_chunk_last=5406,
    )

    assert chunk_trials == pytest.approx(5406)


def test_invalid_chunk_trials():
    with pytest.raises(ValueError, match="no events found inside chunk_range"):
        readers._chunk_attempted_events(
            chunk_range=[10, 20],
            attempted_before_chunk=5,
            attempted_in_chunk_last=None,
        )

    with pytest.raises(ValueError, match="non-positive attempted-event count"):
        readers._chunk_attempted_events(
            chunk_range=[10, 20],
            attempted_before_chunk=100,
            attempted_in_chunk_last=100,
        )


# Check that exactly constant first weights do not trigger the weighted-sample flag
def test_weight_summary_accepts_constant_weights(hepmc_file):
    summary = readers.read_hepmc3_weight_summary(hepmc_file([2.5, 2.5, 2.5]))

    assert summary["nevents"] == 3
    assert not summary["is_nonconstant"]
    assert summary["reference_weight"] == pytest.approx(2.5)
    assert summary["min_weight"] == pytest.approx(2.5)
    assert summary["max_weight"] == pytest.approx(2.5)
    assert summary["finite_weight_count"] == 3
    assert summary["is_constant"]
    assert summary["last_xsection_pb"] == pytest.approx(3.0)
    assert summary["last_xsection_pb_err"] == pytest.approx(0.3)


# Check that non-constant first weights are detected with event-level diagnostics
def test_weight_summary_detects_nonconstant_weights(hepmc_file):
    summary = readers.read_hepmc3_weight_summary(hepmc_file([2.5, 2.5, 2.5000001]))

    assert summary["nevents"] == 3
    assert summary["is_nonconstant"]
    assert not summary["is_constant"]
    assert summary["reference_weight"] == pytest.approx(2.5)
    assert summary["first_nonconstant_event"] == 2
    assert summary["first_nonconstant_weight"] == pytest.approx(2.5000001)
    assert summary["max_abs_deviation"] == pytest.approx(1.0e-7)


# Check that the scan honors the same maxevents boundary used by iceplot
def test_weight_summary_honors_maxevents(hepmc_file):
    summary = readers.read_hepmc3_weight_summary(hepmc_file([1.0, 3.0, 5.0]), maxevents=1)

    assert summary["nevents"] == 1
    assert summary["maxevents_reached"]
    assert not summary["is_nonconstant"]
    assert summary["is_constant"]
    assert summary["reference_weight"] == pytest.approx(1.0)


# Check the full-reader fallback retains the terminal cross section at its scan limit
def test_weight_summary_reader_terminal_xsection(hepmc_file):
    xsections = [(10.0, 1.0), (20.0, 2.0), (30.0, 3.0)]
    path = hepmc_file(xsections=xsections)
    full_summary = readers.read_hepmc3_weight_summary(path)
    limited_summary = readers.read_hepmc3_weight_summary(path, maxevents=2)

    assert full_summary["last_xsection_pb"] == pytest.approx(30.0)
    assert full_summary["last_xsection_pb_err"] == pytest.approx(3.0)
    assert limited_summary["last_xsection_pb"] == pytest.approx(20.0)
    assert limited_summary["last_xsection_pb_err"] == pytest.approx(2.0)


# Check the official reader retains the full weight summary
@pytest.mark.parametrize("compressed", [False, True])
def test_read_hepmc3_weight_summary(tmp_path, compressed):
    path = tmp_path / "weights.hepmc3"
    write_muon_events(path, [2.5, 2.5, 2.5000001])
    if compressed:
        packed = path.with_suffix('.hepmc3.gz')
        packed.write_bytes(gzip.compress(path.read_bytes()))
        path = packed
    summary = readers.read_hepmc3_weight_summary(str(path))

    assert summary["nevents"] == 3
    assert summary["finite_weight_count"] == 3
    assert summary["reference_weight"] == pytest.approx(2.5)
    assert summary["min_weight"] == pytest.approx(2.5)
    assert summary["max_weight"] == pytest.approx(2.5000001)
    assert summary["mean_weight"] == pytest.approx((2.5 + 2.5 + 2.5000001) / 3.0)
    assert summary["first_nonconstant_event"] == 2
    assert summary["is_nonconstant"]


# Check maximum weight overflow attributes select header normalization
def test_weight_summary_detects_overflow_ratios(tmp_path):
    path = tmp_path / "overflow.hepmc3"
    write_muon_events(
        path,
        [1.0, 2.5, 1.0],
        xsections=[(10.0, 1.0), (10.0, 1.0), (10.0, 1.0)],
        overflow_indices={1},
    )

    summary = readers.read_hepmc3_weight_summary(str(path))
    resolved, reason = readers.resolve_hepmc3_xsmode("auto", summary)

    assert summary["maximum_weight_overflow_count"] == 1
    assert summary["is_nonconstant"]
    assert resolved == "header"
    assert "maximum weight overflow ratios" in reason


# Reject an incomplete event instead of accepting its weight and cross section
def test_weight_summary_rejects_incomplete_event(tmp_path):
    path = tmp_path / "incomplete.hepmc3"
    write_muon_events(path, [1.0])
    lines = path.read_text().splitlines()
    last_particle = max(index for index, line in enumerate(lines) if line.startswith('P '))
    path.write_text('\n'.join(lines[:last_particle] + lines[last_particle + 1:]) + '\n')

    with pytest.raises(ValueError, match='Incomplete particle'):
        readers.read_hepmc3_weight_summary(str(path))


# Check header-mode chunks skip directly to their first selected event
def test_hepmc_header_chunk_reader_skip(hepmc_file, reader_calls):
    path = hepmc_file(
        weights=[1.0] * 5, attempted=[2, 5, 9, 14, 20],
        xsections=[(xs, xs / 10.0) for xs in [5.0, 10.0, 15.0, 20.0, 31.0]],
    )

    first_output = readers.read_hepmc3(
        hepmc3file=path,
        obs=[{"event": {"func": lambda event: event.evt.event_number()}}],
        pid=[[13, -13]],
        cuts=[None],
        chunk_range=[0, 1],
        xsmode="header",
        header_xsection=(31.0, 3.1),
    )
    second_output = readers.read_hepmc3(
        hepmc3file=path,
        obs=[{"event": {"func": lambda event: event.evt.event_number()}}],
        pid=[[13, -13]],
        cuts=[None],
        chunk_range=[2, 3],
        xsmode="header",
        header_xsection=(31.0, 3.1),
    )

    assert reader_calls["skip"] == [1]
    assert reader_calls["read"] == [0, 1, 0, 2, 3]
    assert first_output[0]["data"]["event"].tolist() == [0, 1]
    assert second_output[0]["data"]["event"].tolist() == [2, 3]
    assert first_output[0]["xsection_pb"] == pytest.approx(31.0)
    assert second_output[0]["xsection_pb"] == pytest.approx(31.0)
    assert first_output[0]["xsection_pb_err"] == pytest.approx(3.1)
    assert second_output[0]["xsection_pb_err"] == pytest.approx(3.1)


# Check sample-mode chunks retain the preceding attempted-event counter
def test_hepmc_sample_chunk_reader_skip(hepmc_file, reader_calls):
    path = hepmc_file(weights=[1.0] * 5, attempted=[2, 5, 9, 14, 20])

    output = readers.read_hepmc3(
        hepmc3file=path,
        obs=[{"event": {"func": lambda event: event.evt.event_number()}}],
        pid=[[13, -13]],
        cuts=[None],
        chunk_range=[2, 4],
        xsmode="sample",
        header_xsection=(999.0, 99.0),
    )

    assert reader_calls["skip"] == []
    assert reader_calls["read"] == [0, 1, 2, 3, 4]
    assert output[0]["data"]["event"].tolist() == [2, 3, 4]
    assert output[0]["xsection_pb"] == pytest.approx(2.0e11)


# Check chunk partitions preserve physical event weights and cross-section normalization
@pytest.mark.parametrize("xsmode", ["header", "sample"])
@pytest.mark.parametrize("weight_name", [None, "Weight"])
@pytest.mark.parametrize("compressed", [False, True])
def test_hepmc_partition_events_and_xs(hepmc_file, xsmode, weight_name, compressed):
    weights = [1.0, 2.0, 3.0, 4.0, 5.0]
    path = hepmc_file(weights=weights, attempted=[2, 5, 9, 14, 20], weight_name=weight_name)
    if compressed:
        from pathlib import Path
        packed = Path(path + '.gz')
        packed.write_bytes(gzip.compress(Path(path).read_bytes()))
        path = str(packed)
    results = [readers.read_hepmc3(
        hepmc3file=path, obs=[{"event": {"func": lambda event: event.evt.event_number()}}],
        pid=[[13, -13]], cuts=[None], chunk_range=chunk, xsmode=xsmode,
        header_xsection=(31.0, 3.1),
    )[0] for chunk in ([0, 4], [0, 0], [1, 2], [3, 4])]
    full, *parts = results
    np.testing.assert_array_equal(np.concatenate([p["data"]["event"] for p in parts]), np.arange(5))
    np.testing.assert_allclose(np.concatenate([p["weights"] for p in parts]), weights)
    counts = [p["sample_attempted_events"] for p in parts]
    assert sum(counts) == full["sample_attempted_events"] == (20 if xsmode == "sample" else 5)
    sigma = np.average([p["xsection_pb"] for p in parts], weights=counts)
    assert sigma == pytest.approx(full["xsection_pb"])
    assert sigma == pytest.approx(sum(weights) / 20 * 1e12 if xsmode == "sample" else 31.0)


# Check sample-mode uncertainty retains tiny representable event weights
def test_hepmc_error_underflow(hepmc_file):
    path = hepmc_file(weights=[1.0e-200, 2.0e-200])

    output = readers.read_hepmc3(
        hepmc3file=path,
        obs=[{"event": {"func": lambda event: event.evt.event_number()}}],
        pid=[[13, -13]],
        cuts=[None],
        chunk_range=[0, 1],
        xsmode="sample",
    )

    assert output[0]["xsection_pb"] == pytest.approx(1.5e-188, abs=0.0)
    assert output[0]["xsection_pb_err"] == pytest.approx(np.sqrt(5.0) * 0.5e-188, abs=0.0)
    assert output[0]["xsection_pb_err"] > 0.0


# Include rejected attempts in the fiducial MC integral variance
def test_sample_xs_error_includes_rejected_events(hepmc_file, cut_file):
    path = hepmc_file(weights=[1.0e-12, 2.0e-12, 3.0e-12], attempted=[2, 4, 6])
    cuts = cut_file("selected", {}, "return event.evt.event_number() > 0")
    result = readers.read_hepmc3(
        hepmc3file=path, obs=[{}], pid=[[13, -13]], cuts=[cuts], xsmode="sample",
    )[0]
    assert result["xsection_pb"] == pytest.approx(5.0 / 6.0)
    assert result["xsection_pb_err"] == pytest.approx(np.sqrt(13.0) / 6.0)


# Keep zero-weight event samples finite after fiducial selection
def test_zero_weight_sample_is_finite(hepmc_file):
    path = hepmc_file(weights=[0.0, 0.0])
    result = readers.read_hepmc3(
        hepmc3file=path, obs=[{}], pid=[[13, -13]], cuts=[None], xsmode="sample",
    )[0]
    assert result["acceptance"] == pytest.approx(0.0)
    assert result["xsection_pb"] == pytest.approx(0.0)
    assert result["xsection_pb_err"] == pytest.approx(0.0)


# Check the final event remains readable after positioning the last chunk
def test_hepmc_final_single_event_chunk_reader_skip(hepmc_file, reader_calls):
    path = hepmc_file(weights=[1.0] * 5, attempted=[2, 5, 9, 14, 20])

    output = readers.read_hepmc3(
        hepmc3file=path,
        obs=[{"event": {"func": lambda event: event.evt.event_number()}}],
        pid=[[13, -13]],
        cuts=[None],
        chunk_range=[4, 4],
        xsmode="header",
    )

    assert reader_calls["skip"] == [3]
    assert reader_calls["read"] == [0, 4]
    assert output[0]["data"]["event"].tolist() == [4]


# Check automatic reading keeps overflow ratios and uses the header cross section
def test_hepmc_auto_overflow(hepmc_file):
    path = hepmc_file(
        weights=[1.0, 2.5, 1.0], xsections=[(12.0, 0.6)] * 3, overflow_indices={1},
    )

    output = readers.read_hepmc3(
        hepmc3file=path,
        obs=[{"event": {"func": lambda event: event.evt.event_number()}}],
        pid=[[13, -13]],
        cuts=[None],
        chunk_range=[0, 2],
        xsmode="auto",
    )

    assert output[0]["weights"].tolist() == pytest.approx([1.0, 2.5, 1.0])
    assert output[0]["sample_xsection_pb"] == pytest.approx(12.0)
    assert output[0]["sample_xsection_pb_err"] == pytest.approx(0.6)


# Check that auto x-section mode uses header mode for finite constant event weights
def test_auto_xs_constant_weights():
    module = load_iceplot_module()
    resolved, reason = module.resolve_xsmode(
        requested_xsmode="auto",
        weight_summary={
            "is_constant": True,
            "is_nonconstant": False,
            "missing_weight_count": 0,
            "nonfinite_weight_count": 0,
        },
    )

    assert resolved == "header"
    assert "constant event weights" in reason


# Check terminal header extraction accepts finite values and rejects invalid data
def test_iceplot_terminal_header_xsection_validation():
    assert readers.terminal_header_xsection(
        {
            "last_xsection_pb": 12.5,
            "last_xsection_pb_err": 0.4,
        }
    ) == pytest.approx((12.5, 0.4))

    with pytest.raises(ValueError, match="last scanned event"):
        readers.terminal_header_xsection({})
    with pytest.raises(ValueError, match="finite cross-section"):
        readers.terminal_header_xsection(
            {
                "last_xsection_pb": np.nan,
                "last_xsection_pb_err": 0.4,
            }
        )


# Check that auto x-section mode uses sample mode for non-constant event weights
def test_auto_xs_variable_weights():
    module = load_iceplot_module()
    resolved, reason = module.resolve_xsmode(
        requested_xsmode="auto",
        weight_summary={
            "is_constant": False,
            "is_nonconstant": True,
            "missing_weight_count": 0,
            "nonfinite_weight_count": 0,
        },
    )

    assert resolved == "sample"
    assert "non-constant event weights" in reason


# Check maximum weight overflow ratios override generic variable weight detection
def test_auto_xs_overflow():
    module = load_iceplot_module()
    resolved, reason = module.resolve_xsmode(
        requested_xsmode="auto",
        weight_summary={
            "is_constant": False,
            "is_nonconstant": True,
            "missing_weight_count": 0,
            "nonfinite_weight_count": 0,
            "maximum_weight_overflow_count": 1,
        },
    )

    assert resolved == "header"
    assert "maximum weight overflow ratios" in reason


# Check that explicit x-section modes remain explicit after weight scanning
def test_iceplot_resolve_xsmode_explicit_modes():
    module = load_iceplot_module()

    assert module.resolve_xsmode("header", {"is_nonconstant": True})[0] == "header"
    assert module.resolve_xsmode("sample", {"is_constant": True})[0] == "sample"


# Check that auto mode refuses unusable event-weight scans
def test_auto_xs_invalid_weights():
    module = load_iceplot_module()

    with pytest.raises(ValueError, match="no first weight"):
        module.resolve_xsmode(
            requested_xsmode="auto",
            weight_summary={"missing_weight_count": 1, "nonfinite_weight_count": 0},
        )

    with pytest.raises(ValueError, match="non-finite"):
        module.resolve_xsmode(
            requested_xsmode="auto",
            weight_summary={"missing_weight_count": 0, "nonfinite_weight_count": 1},
        )


def test_hepmc_skips_per_event_reco_failures(capsys, hepmc_file):
    output, messages = run_reconstruction_failure_case(capsys, hepmc_file, verbose=1)

    assert output[0]["data"]["stable"].tolist() == [2]
    assert output[0]["data"]["flaky"].tolist() == [12]
    assert output[0]["weights"].tolist() == [1.0]
    assert output[0]["acceptance"] == pytest.approx(1.0 / 3.0)
    assert all("Cannot reconstruct observable" not in message for message in messages)
    assert not any("observable reconstruction failures" in message for message in messages)


def test_hepmc_reco_failure_count(capsys, hepmc_file):
    output, messages = run_reconstruction_failure_case(capsys, hepmc_file, verbose=2)

    assert output[0]["data"]["stable"].tolist() == [2]
    assert output[0]["data"]["flaky"].tolist() == [12]
    assert output[0]["acceptance"] == pytest.approx(1.0 / 3.0)
    assert all("Cannot reconstruct observable" not in message for message in messages)
    assert any("observable reconstruction failures" in message for message in messages)
    assert any('observable="flaky"' in message and "skipped=2" in message for message in messages)
    assert any(
        "first_event=0" in message and "missing cascade for event 0" in message
        for message in messages
    )
