# HepMC3 and HEPData input/output functions
#
# Run tests with: pytest tests/technical/data/test_readers_utils.py -s
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import os
from collections.abc import Generator
from types import ModuleType

import numpy as np
from pyHepMC3 import HepMC3 as hepmc3
from termcolor import cprint

from core.io import steering
from core.io.cache import EventCache, freeze
from core.io.hepdata_reader import read as read_config_hepdata
from core.plot import plot
from core.stats.hist import bins2binwidth, hist
from core.stats.uncertainty import sample_cross_section, stable_sum, stable_sum_squares


# Object wrapper (for cache)
class ice_event:
    # Initialize one selection view of a physical event
    def __init__(self, evt, pid, cut_param, event_cache=None):
        self.evt = evt
        self.pid = pid
        self.cut_param = cut_param
        self._icecache = EventCache() if event_cache is None else event_cache


# Compute one exact cache signature for a particle and cut-parameter view
def _selection_cache_signature(pid, selection: ModuleType):
    return freeze(pid, strict=True), freeze(selection.cut_param, strict=True)


# Assign cache ids to selections with equivalent particle and cut parameters
def _selection_cache_ids(pid: list, selections: list[ModuleType]) -> tuple[list[int], int, list]:
    if len(pid) != len(selections):
        raise ValueError("readers: PID and selection counts must agree")

    ids = []
    selection_signatures = []
    signatures = {}
    next_id = 0
    for selection_id in range(len(selections)):
        try:
            signature = _selection_cache_signature(pid=pid[selection_id], selection=selections[selection_id])
        except (RecursionError, TypeError, ValueError):
            ids.append(next_id)
            selection_signatures.append(None)
            next_id += 1
            continue

        cache_id = signatures.get(signature)
        if cache_id is None:
            cache_id = next_id
            signatures[signature] = cache_id
            next_id += 1
        ids.append(cache_id)
        selection_signatures.append(signature)

    counts = [0] * next_id
    for cache_id in ids:
        counts[cache_id] += 1
    shared_signatures = [
        signature if signature is not None and counts[cache_id] > 1 else None
        for cache_id, signature in zip(ids, selection_signatures, strict=True)
    ]
    return ids, next_id, shared_signatures


# Require shared particle and cut parameters to remain immutable during evaluation
def _validate_selection_cache(pid, selection: ModuleType, expected) -> None:
    if expected is None:
        return
    try:
        current = _selection_cache_signature(pid=pid, selection=selection)
    except (RecursionError, TypeError, ValueError) as exc:
        raise RuntimeError("readers: A shared event selection changed to unsupported PID or cut parameters") from exc
    if current != expected:
        raise RuntimeError("readers: Shared event PID and cut parameters must not change during evaluation")


# Load and validate a built-in or file-based event-selection module
def load_cut_module(reference: str | None) -> ModuleType:
    resolved = "core.analysis.cuts.default" if reference is None else reference
    module = steering.load_python_module(resolved)

    if not callable(getattr(module, "cut_func", None)):
        raise AttributeError(f"Cut definition '{reference}' must define callable cut_func(event)")
    if not hasattr(module, "cut_param"):
        raise AttributeError(f"Cut definition '{reference}' must define cut_param")
    return module


# Construct every requested observable for one accepted event view
def _project_observables(event, observables: dict, failures: dict, event_id: int, verbosity: int):
    values = {}
    for name, observable in observables.items():
        try:
            values[name] = observable["func"](event)
        except Exception as exc:
            failure = failures[name]
            failure["count"] += 1
            if failure["first_event"] is None:
                failure["first_event"] = event_id
                failure["first_error"] = str(exc)
            if verbosity >= 3 and len(failure["examples"]) < 3:
                failure["examples"].append((event_id, str(exc)))
            return None
    return values


def clean_filename(txt: str) -> str:
    """Clean filename for saving"""
    txt = txt.replace(" ]", "]").replace("[ ", "[").replace("  ", " ").replace(" ", "-")
    txt = txt.replace("\\;", "").replace("\\", "").replace("$", "").replace("^", "")
    txt = txt.replace("{", "(").replace("}", ")")

    return txt


def get_observables(module_name: str) -> dict:
    """Read out observables constructions module under 'core.analysis.observables'"""
    reference = steering.resolve_python_reference(
        module_name, package="core.analysis.observables", dataset_path=None, cdir=os.getcwd()
    )
    module = steering.load_python_module(reference)
    return steering.observables_from_module(module, reference=reference)


def parallel_wrapper(param: dict) -> dict:
    """Wrapper for parallel (multiprocess) processing"""

    ## Read event data
    mcdata = read_hepmc3(
        chunk_range=param["chunk_range"],
        hepmc3file=param["hepmc3file"],
        obs=param["obs"],
        pid=param["pid"],
        cuts=param["cuts"],
        xsmode=param["xsmode"],
        header_xsection=param["header_xsection"],
        verbose=param["verbose"],
    )

    ## Histograms per dataset
    h = [None] * len(mcdata)
    scales = param["scales"]
    if len(scales) != len(h):
        raise ValueError("parallel_wrapper: scales must provide one value per dataset set")

    for i in range(len(h)):
        h[i] = plot.histmc(
            mcdata=mcdata[i],
            obs=param["obs"][i],
            scale=scales[i],
            density=param["density"],
            density_uncertainty=param["density_uncertainty"],
            covariance_mode=param.get("covariance_mode", "diagonal"),
            color=plot.colors(param["k"]),
            label=param["label"],
        )

    if not param.get("return_diagnostics", False):
        return h
    return {"histograms": h, "diagnostics": hepmc3_chunk_diagnostics(mcdata=mcdata, filename=param["hepmc3file"])}


# Compute additive before-cut and after-cut statistics for one HepMC3 chunk
def hepmc3_chunk_diagnostics(mcdata: list[dict], filename: str) -> dict:
    sample_weights = np.asarray(mcdata[0]["sample_weights"], dtype=float)
    return {
        "filename": filename,
        "before": {
            "nevents": len(sample_weights),
            "attempted_events": mcdata[0]["sample_attempted_events"],
            "wsum": stable_sum(sample_weights),
            "wsum2": stable_sum_squares(sample_weights),
            "xsection_pb": mcdata[0]["sample_xsection_pb"],
            "xsection_pb_err": mcdata[0]["sample_xsection_pb_err"],
        },
        "after": [
            {
                "nevents": len(record["weights"]),
                "wsum": stable_sum(record["weights"]),
                "wsum2": stable_sum_squares(record["weights"]),
            }
            for record in mcdata
        ],
    }


def _chunk_attempted_events(
    *, chunk_range: list, attempted_before_chunk: float, attempted_in_chunk_last: float
) -> float:
    """Return the exact number of attempted events represented by a chunk"""
    if attempted_in_chunk_last is None:
        raise ValueError(f"_chunk_attempted_events: no events found inside chunk_range = {chunk_range}")

    chunk_trials = attempted_in_chunk_last - attempted_before_chunk
    if chunk_trials <= 0:
        raise ValueError(
            f"_chunk_attempted_events: non-positive attempted-event count {chunk_trials} for chunk_range = {chunk_range}"
        )
    return chunk_trials


def _verbosity_level(verbose) -> int:
    try:
        return int(verbose)
    except (TypeError, ValueError):
        return 0


# Require positive total, chunk and minimum partition sizes
def _validate_partitions(name, total, chunk, minimum, automatic=False) -> None:
    labels = ("Maximum chunks" if automatic else "Chunk size", "Minimum chunk size", "Total size")
    for label, value in zip(labels, (chunk, minimum, total), strict=True):
        if value <= 0:
            raise ValueError(f"{name}: {label} must be greater than 0.")


# Generate fixed-size inclusive ranges while controlling the trailing chunk
def generate_partitions(
    totalsize: int, chunksize: int, min_chunksize: int = 1
) -> Generator[tuple[int, int], None, None]:
    _validate_partitions("generate_partitions", totalsize, chunksize, min_chunksize)

    for start in range(0, totalsize, chunksize):
        stop = min(start + chunksize, totalsize)
        if totalsize - stop < min_chunksize:
            stop = totalsize
        yield start, stop - 1
        if stop == totalsize:
            break


# Generate balanced CPU-aware ranges without sub-minimum chunks
def generate_auto_partitions(
    totalsize: int, max_chunks: int, min_chunksize: int = 1000
) -> Generator[tuple[int, int], None, None]:
    _validate_partitions("generate_auto_partitions", totalsize, max_chunks, min_chunksize, True)

    chunk_count = min(max_chunks, max(1, totalsize // min_chunksize))
    base_size, extra_events = divmod(totalsize, chunk_count)

    start = 0
    for index in range(chunk_count):
        size = base_size + (1 if index < extra_events else 0)
        end = start + size - 1
        yield (start, end)
        start = end + 1


def generate_chunks(lst: list, n: int) -> Generator[list, None, None]:
    """Yield successive n-sized chunks from the list"""
    if n <= 0:
        raise ValueError("generate_chunks: Chunk size must be greater than zero.")

    for i in range(0, len(lst), n):
        yield lst[i : i + n]


# Compute an empty HepMC3 event-count and weight summary
def _new_hepmc3_weight_summary(rel_tol: float, abs_tol: float) -> dict:
    return {
        "nevents": 0,
        "maxevents_reached": False,
        "reference_weight": None,
        "first_nonconstant_event": None,
        "first_nonconstant_weight": None,
        "min_weight": None,
        "max_weight": None,
        "mean_weight": None,
        "max_abs_deviation": 0.0,
        "finite_weight_count": 0,
        "missing_weight_count": 0,
        "nonfinite_weight_count": 0,
        "maximum_weight_overflow_count": 0,
        "lhe_weight_count": 0,
        "last_xsection_pb": None,
        "last_xsection_pb_err": None,
        "is_constant": False,
        "is_nonconstant": False,
        "rel_tol": rel_tol,
        "abs_tol": abs_tol,
    }


# Compute whether two finite weights differ beyond the scan tolerances
def _weights_differ(first: float, second: float, rel_tol: float, abs_tol: float) -> bool:
    scale = max(abs(first), abs(second))
    return abs(first - second) > abs_tol + rel_tol * scale


# Combine per-file weight scans into one logical HepMC3 source summary
def combine_hepmc3_weight_summaries(summaries: list[dict]) -> dict:
    if not summaries:
        raise ValueError("readers: cannot combine an empty list of HepMC3 summaries")

    rel_tol = summaries[0]["rel_tol"]
    abs_tol = summaries[0]["abs_tol"]
    combined = _new_hepmc3_weight_summary(rel_tol=rel_tol, abs_tol=abs_tol)
    combined["_weight_sum"] = 0.0

    for summary in summaries:
        same_rel_tol = np.isclose(summary["rel_tol"], rel_tol, rtol=0.0, atol=0.0)
        same_abs_tol = np.isclose(summary["abs_tol"], abs_tol, rtol=0.0, atol=0.0)
        if not same_rel_tol or not same_abs_tol:
            raise ValueError("readers: HepMC3 weight summaries use different tolerances")

        event_offset = combined["nevents"]
        reference = summary.get("reference_weight")
        if combined["reference_weight"] is None and reference is not None:
            combined["reference_weight"] = reference

        global_reference = combined["reference_weight"]
        local_first = summary.get("first_nonconstant_event")
        if combined["first_nonconstant_event"] is None:
            if local_first is not None:
                combined["first_nonconstant_event"] = event_offset + local_first
                combined["first_nonconstant_weight"] = summary.get("first_nonconstant_weight")
            elif (
                reference is not None
                and global_reference is not None
                and _weights_differ(reference, global_reference, rel_tol, abs_tol)
            ):
                combined["first_nonconstant_event"] = event_offset
                combined["first_nonconstant_weight"] = reference

        for key in ("nevents", "finite_weight_count", "missing_weight_count", "nonfinite_weight_count"):
            combined[key] += summary[key]
        combined["maximum_weight_overflow_count"] += summary.get("maximum_weight_overflow_count", 0)
        combined["lhe_weight_count"] += summary.get("lhe_weight_count", 0)
        combined["maxevents_reached"] |= summary["maxevents_reached"]
        combined["_weight_sum"] += float(summary.get("mean_weight") or 0.0) * summary["finite_weight_count"]

        for key, reducer in (("min_weight", min), ("max_weight", max)):
            value = summary.get(key)
            if value is not None:
                current = combined[key]
                combined[key] = value if current is None else reducer(current, value)

        if global_reference is not None:
            deviations = [
                abs(value - global_reference)
                for value in (summary.get("min_weight"), summary.get("max_weight"))
                if value is not None
            ]
            if deviations:
                combined["max_abs_deviation"] = max(combined["max_abs_deviation"], max(deviations))

        if summary["nevents"] > 0:
            combined["last_xsection_pb"] = summary.get("last_xsection_pb")
            combined["last_xsection_pb_err"] = summary.get("last_xsection_pb_err")

    return _finalize_hepmc3_weight_summary(combined)


# Add one event weight to a HepMC3 summary
def _accumulate_hepmc3_event_weight(summary: dict, weight: float | None) -> None:
    event_index = summary["nevents"]

    if weight is None:
        summary["missing_weight_count"] += 1
        if summary["reference_weight"] is not None and summary["first_nonconstant_event"] is None:
            summary["first_nonconstant_event"] = event_index
        summary["nevents"] += 1
        return

    if not np.isfinite(weight):
        summary["nonfinite_weight_count"] += 1
        if summary["reference_weight"] is not None and summary["first_nonconstant_event"] is None:
            summary["first_nonconstant_event"] = event_index
            summary["first_nonconstant_weight"] = weight
        summary["nevents"] += 1
        return

    if summary["reference_weight"] is None:
        summary["reference_weight"] = weight

    summary["finite_weight_count"] += 1
    summary["_weight_sum"] += weight
    summary["min_weight"] = weight if summary["min_weight"] is None else min(summary["min_weight"], weight)
    summary["max_weight"] = weight if summary["max_weight"] is None else max(summary["max_weight"], weight)

    abs_deviation = abs(weight - summary["reference_weight"])
    summary["max_abs_deviation"] = max(summary["max_abs_deviation"], abs_deviation)

    scale = max(abs(weight), abs(summary["reference_weight"]))
    tolerance = summary["abs_tol"] + summary["rel_tol"] * scale
    if abs_deviation > tolerance and summary["first_nonconstant_event"] is None:
        summary["first_nonconstant_event"] = event_index
        summary["first_nonconstant_weight"] = weight

    summary["nevents"] += 1


# Retain the cross section carried by the last complete scanned event
def _record_hepmc3_event_cross_section(summary: dict, xsection_pb: float | None, xsection_pb_err: float | None) -> None:
    summary["last_xsection_pb"] = xsection_pb
    summary["last_xsection_pb_err"] = xsection_pb_err


# Add one complete event's weight, cross section and overflow marker to a HepMC3 summary
def _record_hepmc3_event(
    summary: dict,
    weight: float | None,
    xsection_pb: float | None,
    xsection_pb_err: float | None,
    maximum_weight_overflow: bool,
) -> None:
    _accumulate_hepmc3_event_weight(summary, weight)
    _record_hepmc3_event_cross_section(summary, xsection_pb, xsection_pb_err)
    summary["maximum_weight_overflow_count"] += int(maximum_weight_overflow)


# Compute whether one HepMC3 event carries a named attribute
def _hepmc3_event_has_attribute(evt, name: str) -> bool:
    attribute_names = getattr(evt, "attribute_names", None)
    if not callable(attribute_names):
        return False
    return name in list(attribute_names())


# Validate one shared terminal HepMC3 header cross section
def validate_header_xsection(header_xsection) -> tuple[float, float]:
    if not isinstance(header_xsection, (tuple, list)) or len(header_xsection) != 2:
        raise ValueError("validate_header_xsection: expected a (cross section, uncertainty) pair")

    xsection_pb = float(header_xsection[0])
    xsection_pb_err = float(header_xsection[1])
    if not np.isfinite(xsection_pb) or not np.isfinite(xsection_pb_err):
        raise ValueError("validate_header_xsection: expected finite cross-section values")
    return xsection_pb, xsection_pb_err


# Compute the terminal finite HepMC3 header cross section found by a pre-scan
def terminal_header_xsection(weight_summary: dict) -> tuple[float, float]:
    xsection_pb = weight_summary.get("last_xsection_pb")
    xsection_pb_err = weight_summary.get("last_xsection_pb_err")
    if xsection_pb is None or xsection_pb_err is None:
        raise ValueError("terminal_header_xsection: GenCrossSection is missing from the last scanned event")
    return validate_header_xsection((xsection_pb, xsection_pb_err))


# Complete the derived fields of a HepMC3 weight summary
def _finalize_hepmc3_weight_summary(summary: dict) -> dict:
    if summary["finite_weight_count"] > 0:
        summary["mean_weight"] = summary["_weight_sum"] / summary["finite_weight_count"]

    summary["is_nonconstant"] = (
        summary["first_nonconstant_event"] is not None
        or summary["missing_weight_count"] > 0
        or summary["nonfinite_weight_count"] > 0
    )
    summary["is_constant"] = (
        summary["nevents"] > 0
        and summary["finite_weight_count"] == summary["nevents"]
        and not summary["is_nonconstant"]
    )
    summary.pop("_weight_sum", None)
    return summary


# Read event weights and cross sections through the official HepMC3 reader
def read_hepmc3_weight_summary(
    filename: str, maxevents: int = int(1e9), rel_tol: float = 1e-12, abs_tol: float = 0.0
) -> dict:
    summary = _new_hepmc3_weight_summary(rel_tol=rel_tol, abs_tol=abs_tol)
    summary["_weight_sum"] = 0.0
    file = hepmc3.deduce_reader(filename)
    try:
        while not file.failed():
            if summary["nevents"] == maxevents:
                summary["maxevents_reached"] = True
                break

            evt = hepmc3.GenEvent()
            file.read_event(evt)

            if file.failed():  # end-of-file
                break

            # A missing particle record can leave an uninitialized GenParticle in HepMC3
            if any(p.pid() == 0 and p.status() == 0 for p in evt.particles()):
                raise ValueError(f'readers: Incomplete particle in {filename}, event {evt.event_number()}')

            weights = evt.weights()
            cross_section = evt.cross_section()
            summary["lhe_weight_count"] += int(_hepmc3_event_has_attribute(evt, "graniitti_lhe_weight"))
            _record_hepmc3_event(
                summary,
                None if len(weights) == 0 else float(weights[0]),
                None if cross_section is None else float(cross_section.xsec()),
                None if cross_section is None else float(cross_section.xsec_err()),
                _hepmc3_event_has_attribute(evt, "maximum_weight_overflow"),
            )

    finally:
        file.close()
    return _finalize_hepmc3_weight_summary(summary)


# Resolve automatic HepMC3 cross-section normalization from the event weights
def resolve_hepmc3_xsmode(requested_xsmode: str, weight_summary: dict) -> tuple[str, str]:
    if requested_xsmode in ("header", "sample"):
        return requested_xsmode, f"user requested {requested_xsmode}"
    if requested_xsmode != "auto":
        raise ValueError(f"readers: Unknown xsmode = {requested_xsmode}")
    if weight_summary.get("missing_weight_count", 0) > 0:
        raise ValueError("readers: Cannot resolve xsmode='auto' because at least one scanned event has no first weight")
    if weight_summary.get("nonfinite_weight_count", 0) > 0:
        raise ValueError("readers: Cannot resolve xsmode='auto' because at least one scanned first weight is non-finite")
    if weight_summary.get("lhe_weight_count", 0) > 0:
        if weight_summary["lhe_weight_count"] != weight_summary["nevents"]:
            raise ValueError("readers: Cannot combine native and LHE event weight conventions")
        return "header", "auto selected header from LHE event weights"
    if weight_summary.get("maximum_weight_overflow_count", 0) > 0:
        return "header", "auto selected header from maximum weight overflow ratios"
    if weight_summary.get("is_nonconstant", False):
        return "sample", "auto selected sample from non-constant event weights"
    if weight_summary.get("is_constant", False):
        return "header", "auto selected header from constant event weights"
    raise ValueError("readers: Cannot resolve xsmode='auto' from an empty or invalid weight scan")


# Skip an exact number of events with ReaderAscii
def _skip_hepmc3_events(file, event_count: int) -> None:
    if event_count < 0:
        raise ValueError("_skip_hepmc3_events: event count must be nonnegative")
    if event_count == 0:
        return

    # Read the run info before skip(), which can otherwise consume the first event header
    file.read_event(hepmc3.GenEvent())
    skip_succeeded = not file.failed() and (event_count == 1 or file.skip(event_count - 1))
    if not skip_succeeded or file.failed():
        raise RuntimeError(f"_skip_hepmc3_events: failed to skip {event_count} events")


# Position a ReaderAscii instance at the first event of one chunk
def _position_hepmc3_reader(file, chunk_start: int, xsmode: str) -> float:
    if chunk_start < 0:
        raise ValueError("_position_hepmc3_reader: chunk start must be nonnegative")
    if chunk_start == 0:
        return 0.0

    event_count = chunk_start - 1 if xsmode == "sample" else chunk_start
    _skip_hepmc3_events(file=file, event_count=event_count)

    if xsmode != "sample":
        return 0.0

    previous_event = hepmc3.GenEvent()
    file.read_event(previous_event)
    if file.failed():
        raise RuntimeError(f"_position_hepmc3_reader: failed to read the event preceding chunk start {chunk_start}")
    return float(previous_event.cross_section().get_attempted_events())


# Read HepMC3 events and project the requested fiducial observables
def read_hepmc3(
    hepmc3file: str,
    obs: list,
    pid: list,
    cuts: list,
    chunk_range=None,
    xsmode: str = "auto",
    header_xsection: tuple[float, float] | None = None,
    MC_WEIGHT_RESCALE: float = 1e12,
    verbose: int | bool = False,
) -> dict:
    """Read in a single HepMC3 file and compute observables"""
    chunk_range = [0, int(1e9)] if chunk_range is None else chunk_range
    verbosity = _verbosity_level(verbose)

    requested_xsmode = xsmode
    if requested_xsmode == "auto":
        maxevents = int(chunk_range[1]) + 1
        weight_summary = read_hepmc3_weight_summary(filename=hepmc3file, maxevents=maxevents)
        xsmode, reason = resolve_hepmc3_xsmode(requested_xsmode=requested_xsmode, weight_summary=weight_summary)
        if verbosity > 0:
            print(f"{__name__}.read_hepmc3: xsmode={xsmode} ({reason})")

    if verbosity > 0:
        print(
            __name__
            + f".read_hepmc3: {hepmc3file} | chunk_range = {chunk_range} | pid = {pid}, cuts = {cuts}, xsmode = {xsmode}"
        )
    # Per triplet (obs,pid,cuts)

    output = [None] * len(obs)
    MySelection = [None] * len(obs)
    reco_failures = [None] * len(obs)

    for k in range(len(obs)):
        meta = {
            "source": os.path.abspath(hepmc3file),
            "data": {},
            "weights": [],
            "event_ids": [],
            "sample_event_ids": [],
            "sample_weights": [],
            "sample_attempted_events": 0,
            "sample_xsection_pb": 0,
            "sample_xsection_pb_err": 0,
            "xsection_pb": 0,
            "xsection_pb_err": 0,
            "acceptance": 0,
            "event_acceptance": 0,
        }
        for OBS in obs[k]:
            meta["data"][OBS] = []

        output[k] = copy.deepcopy(meta)
        reco_failures[k] = {
            OBS: {"count": 0, "first_event": None, "first_error": None, "examples": []} for OBS in obs[k]
        }

        MySelection[k] = load_cut_module(cuts[k])

    cache_ids, cache_count, cache_signatures = _selection_cache_ids(pid=pid, selections=MySelection)
    # Global

    mc_weights_tot = []
    mc_xsec_pb = 0
    mc_xsec_pb_err = 0
    attempted_before_chunk = 0.0
    attempted_in_chunk_last = None
    index = []
    # Open HepMC3 file
    file = hepmc3.deduce_reader(hepmc3file)
    try:
        attempted_before_chunk = _position_hepmc3_reader(file=file, chunk_start=chunk_range[0], xsmode=xsmode)
        k = chunk_range[0] - 1

        # We read event-by-event
        while not file.failed() and k < chunk_range[1]:
            k += 1
            evt = hepmc3.GenEvent()
            file.read_event(evt)

            if file.failed():  # end-of-file
                break

            # Convert HepMC momenta before applying cuts and projecting observables
            evt.set_units(hepmc3.Units.GEV, hepmc3.Units.MM)

            # Keep track of the cumulative attempted-event counter immediately before the chunk
            current_trials = float(evt.cross_section().get_attempted_events())

            index.append(k)

            # Event weight for each event
            w = float(evt.weights()[0])  # Take the first weight with [0]
            mc_weights_tot.append(w)

            # Get GenCrossSection info attributes
            mc_xsec_pb = float(evt.cross_section().xsec())
            mc_xsec_pb_err = float(evt.cross_section().xsec_err())
            attempted_in_chunk_last = current_trials
            ## Now loop over different (obs, pid, cuts) triplets

            # Share cached physics only between equivalent particle and cut-parameter views
            event_caches = [EventCache() for _ in range(cache_count)]

            for j in range(len(output)):
                event = ice_event(
                    evt=evt, pid=pid[j], cut_param=MySelection[j].cut_param, event_cache=event_caches[cache_ids[j]]
                )
                expected_signature = cache_signatures[j]
                try:
                    # Apply the selected analysis cuts
                    if not MySelection[j].cut_func(event):
                        continue
                    values = _project_observables(
                        event=event, observables=obs[j], failures=reco_failures[j], event_id=k, verbosity=verbosity
                    )
                    if values is None:
                        continue

                    # Collect observables for the accepted event
                    for OBS in obs[j]:
                        output[j]["data"][OBS].append(values[OBS])

                    # Event weight for passed events
                    output[j]["weights"].append(w)
                    output[j]["event_ids"].append(k)
                finally:
                    _validate_selection_cache(pid=pid[j], selection=MySelection[j], expected=expected_signature)
    finally:
        file.close()  #! Close manually (no context manager support in pyHepMC3 library)

    if verbosity >= 2:
        for j in range(len(reco_failures)):
            for OBS, failure in reco_failures[j].items():
                if failure["count"] == 0:
                    continue
                cprint(
                    __name__
                    + f".read_hepmc3: observable reconstruction failures | chunk={chunk_range} | set={j} | "
                    + f'observable="{OBS}" | skipped={failure["count"]} | first_event={failure["first_event"]} | '
                    + f"first_error=[{failure['first_error']}]",
                    "yellow",
                )
                if verbosity >= 3:
                    for event_id, error in failure["examples"]:
                        cprint(
                            __name__
                            + f".read_hepmc3: reconstruction failure example | chunk={chunk_range} | set={j} | "
                            + f'observable="{OBS}" | event={event_id} | error=[{error}]',
                            "yellow",
                        )
    mc_weights_tot = np.asarray(mc_weights_tot)

    # Compute cross section from the sample weights (works with weighted events)
    if xsmode == "sample":
        N = _chunk_attempted_events(
            chunk_range=chunk_range,
            attempted_before_chunk=attempted_before_chunk,
            attempted_in_chunk_last=attempted_in_chunk_last,
        )
        mc_xsec_pb, mc_xsec_pb_err = sample_cross_section(mc_weights_tot, N, MC_WEIGHT_RESCALE)

    elif xsmode == "header":
        if header_xsection is not None:
            mc_xsec_pb, mc_xsec_pb_err = validate_header_xsection(header_xsection)
    else:
        raise Exception(__name__ + f".read_hepmc3: Unknown xsmode = {xsmode}")

    if verbosity > 0:
        print(hepmc3file)
        print(
            f"Before python cuts | mc_xs_pb: {mc_xsec_pb:.3g} +- {mc_xsec_pb_err:.3g} | "
            + f"events read: {len(mc_weights_tot)} | wsum: {np.sum(mc_weights_tot):0.3E} | chunk: {chunk_range}"
        )
    ## Now loop over different (obs, pid, cuts) triplets

    wsum_tot = stable_sum(mc_weights_tot)

    for j in range(len(output)):
        # Convert observables to numpy
        for OBS in output[j]["data"]:
            output[j]["data"][OBS] = np.asarray(output[j]["data"][OBS])

        # Convert weights into numpy
        output[j]["weights"] = np.asarray(output[j]["weights"])
        output[j]["event_ids"] = np.asarray(output[j]["event_ids"], dtype=np.int64)
        output[j]["sample_event_ids"] = np.asarray(index, dtype=np.int64)
        output[j]["sample_weights"] = np.asarray(mc_weights_tot, dtype=float)
        output[j]["sample_attempted_events"] = N if xsmode == "sample" else len(mc_weights_tot)
        output[j]["sample_xsection_pb"] = mc_xsec_pb
        output[j]["sample_xsection_pb_err"] = mc_xsec_pb_err
        wsum = stable_sum(output[j]["weights"])

        # Python cuts fiducial acceptance and cross section
        output[j]["acceptance"] = wsum / wsum_tot if abs(wsum_tot) > 0.0 else 0.0
        output[j]["event_acceptance"] = (
            len(output[j]["weights"]) / len(mc_weights_tot) if len(mc_weights_tot) > 0 else 0.0
        )
        if xsmode == "sample":
            output[j]["xsection_pb"], output[j]["xsection_pb_err"] = sample_cross_section(
                output[j]["weights"], N, MC_WEIGHT_RESCALE
            )
        else:
            output[j]["xsection_pb"] = mc_xsec_pb * output[j]["acceptance"]
            output[j]["xsection_pb_err"] = mc_xsec_pb_err * abs(output[j]["acceptance"])

        if verbosity > 0:
            print(
                "After python cuts  | "
                + f"mc_xs_pb: {output[j]['xsection_pb']:.3g} +- {output[j]['xsection_pb_err']:.3g} | "
                + f"events read: {len(output[j]['weights'])} | wsum: {wsum:0.3E} | "
                + f"event acceptance: {output[j]['event_acceptance']:0.3E} | "
                + f"cross section acceptance: {output[j]['acceptance']:0.3E} <{cuts[j]}>"
            )

    return output


def read_hepdata(
    dataset: dict,
    datapath: str,
    datatype: str,
    all_obs: dict,
    cdir: str = None,
    reader: str = None,
    dataset_path: str = None,
) -> tuple[dict, dict]:
    """Read in HepData files"""

    hepdata = {}
    obs = steering.select_histogram_observables(dataset, all_obs)

    # Over different observables
    for histogram in dataset["hist"]:
        filename = steering.resolve_data_reference(
            histogram["file"], datapath=datapath, dataset_path=dataset_path, cdir=cdir
        )
        OBS = histogram["obs"]
        file_filter = histogram.get("file_filter")
        rebin_factor = histogram.get("rebin_factor")

        ### Read out HEPData
        if datatype == "RAW_SCALAR":
            hepdata[OBS] = ReadHEPData_raw(filename=filename, hist_dict=histogram, density=all_obs[OBS]["density"])

        else:
            try:
                hepdata[OBS] = read_config_hepdata(
                    reader,
                    dataset_path=dataset_path,
                    cdir=cdir,
                    filename=filename,
                    file_filter=file_filter,
                    rebin_factor=rebin_factor,
                    hist=histogram,
                    dataset=dataset,
                    all_obs=all_obs,
                    obs=OBS,
                )
            except (AttributeError, FileNotFoundError, ValueError) as exc:
                raise Exception(__name__ + f'.read_hepdata: Invalid bundle reader "{reader}"') from exc

        # Isolate observable metadata from cached reader results
        hepdata[OBS] = copy.deepcopy(hepdata[OBS])

        # y-axis scale
        hepdata[OBS]["scale"] = histogram["scale"]

        # fit importance weight scalar
        hepdata[OBS]["fitw"] = histogram.get("fitw", 1.0)

        # Retain the explicit positive plotting range required by a logarithmic x axis
        if all_obs[OBS].get("xscale", "linear") == "log":
            xlim = np.asarray(all_obs[OBS]["xlim"], dtype=float)
            if xlim.shape != (2,) or not np.all(np.isfinite(xlim)) or xlim[0] <= 0.0 or xlim[1] <= xlim[0]:
                raise ValueError(
                    f'{__name__}.read_hepdata: Observable "{OBS}" requires positive ordered xlim values for logarithmic x scale'
                )
            hepdata[OBS]["xlim"] = copy.deepcopy(all_obs[OBS]["xlim"])

        ### Set the same binning for the observables as HepDATA has --> for MC
        obs[OBS]["bins"] = copy.deepcopy(hepdata[OBS]["bins"])
        obs[OBS]["binwidth"] = copy.deepcopy(hepdata[OBS]["binwidth"])
        obs[OBS]["xlim"] = copy.deepcopy(hepdata[OBS]["xlim"])
        mc_scales = [source["mc_scale"] for source in (hepdata[OBS], histogram) if "mc_scale" in source]
        if len(mc_scales) > 1:
            raise ValueError(f'readers: MC normalization for "{OBS}" is specified by both reader and card')
        if mc_scales:
            mc_bin_scale = np.asarray(mc_scales[0], dtype=float)
            if (
                mc_bin_scale.ndim > 0 and mc_bin_scale.shape != (len(obs[OBS]["bins"]) - 1,)
                or not np.isfinite(mc_bin_scale).all()
                or np.any(mc_bin_scale <= 0.0)
            ):
                raise ValueError(f'readers: Invalid MC bin scales for "{OBS}"')
            obs[OBS]["mc_scale"] = mc_bin_scale.copy()
        if "valid" in hepdata[OBS]:
            obs[OBS]["valid"] = copy.deepcopy(hepdata[OBS]["valid"])

    return hepdata, obs


# Read one unbinned scalar measurement and construct its histogram
def ReadHEPData_raw(filename: str, hist_dict: dict, density: bool) -> dict:
    x = np.loadtxt(filename, dtype=float, ndmin=1)
    bins = np.linspace(hist_dict["xmin"], hist_dict["xmax"], hist_dict["nbins"] + 1)
    counts, errs, bins, cbins = hist(x=x, bins=bins, density=density)

    return {
        "x": cbins,
        "bins": bins,
        "binwidth": bins2binwidth(bins),
        "xlim": [np.min(bins), np.max(bins)],
        "y": counts,
        "y_err": errs,
        "y_err_stat": errs,
        "y_err_syst": np.zeros_like(errs),
        "uncertainties": [
            {
                "name": "sample_statistical",
                "category": "statistical",
                "correlation": "uncorrelated",
                "up": np.asarray(errs, dtype=float),
                "down": np.asarray(errs, dtype=float),
                "provenance": "Raw sample weighted-count uncertainty",
            }
        ],
    }
