#!/usr/bin/env python
#
# iceplot: MC and data differential comparison via HepMC3 & HEPData steering cards
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import ast
import copy
import glob
import json
import math
import multiprocessing
import os
import pathlib
import sys
import time

import matplotlib.pyplot as plt
from termcolor import cprint

from core import __AUTHOR__, __RELEASE__, __version__
from core.io import cli, readers, steering
from core.io.files import ensure_dir
from core.plot import plot as plots
from core.stats import objective, uncertainty
from core.stats.uncertainty import (
    aggregate_header_xsection,
    aggregate_sample_xsection,
    cross_section_ratio,
    histogram_covariance,
    histogram_mc_precision,
)
from core.stats.validation import check_physics, measurement_assessment
from core.tune.runtime import process as iceruntime


# Parse automatic or explicitly sized event chunks
def parse_chunksize(value):
    text = str(value).strip().lower()
    if text == "auto":
        return text

    try:
        chunksize = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--chunksize must be 'auto' or an integer") from exc

    if chunksize < 1000:
        raise argparse.ArgumentTypeError("--chunksize must be at least 1000 events")
    return chunksize


# Verify inclusive ranges cover the requested sample exactly once
def validate_event_partitions(partitions, nevents):
    expected_start = 0
    for start, end in partitions:
        if start != expected_start or end < start:
            raise RuntimeError("validate_event_partitions: non-contiguous or invalid event ranges")
        expected_start = end + 1

    if expected_start != nevents:
        raise RuntimeError(
            f"validate_event_partitions: covered {expected_start} events, expected {nevents}"
        )


# Build complete event ranges and limit workers to the number of ranges
def build_event_partitions(nevents, chunksize, cores):
    if cores <= 0:
        raise ValueError("build_event_partitions: cores must be positive")

    if chunksize == "auto":
        partitions = list(
            readers.generate_auto_partitions(totalsize=nevents, max_chunks=cores, min_chunksize=1000)
        )
    else:
        partitions = list(
            readers.generate_partitions(totalsize=nevents, chunksize=chunksize, min_chunksize=1000)
        )

    validate_event_partitions(partitions=partitions, nevents=nevents)
    return partitions, min(cores, len(partitions))


# Validate optional per-input MC assignments to dataset samples
def validate_mc_hepdata_samples(assignments, dataset_sets, input_count):
    if assignments is None:
        return
    if len(assignments) != input_count:
        raise ValueError("iceplot: --mc-hepdata-sample must provide one name per HepMC3 input")

    samples = {sample for dataset in dataset_sets for sample in steering.set_sample_names(dataset)}
    unknown = set(assignments).difference(samples)
    if unknown:
        raise ValueError(f"iceplot: --mc-hepdata-sample contains unknown names {sorted(unknown)}")


# Compute the MC input indices assigned to one HEPData dataset set
def mc_indices_for_dataset(assignments, dataset, input_count):
    if dataset.get("mc", True) is False:
        return []
    if assignments is None:
        return range(input_count)

    samples = set(steering.set_sample_names(dataset))
    indices = [index for index, assignment in enumerate(assignments) if assignment in samples]
    if input_count > 0 and not indices:
        raise ValueError(f"iceplot: no HepMC3 inputs assigned to dataset samples {sorted(samples)}")
    return indices


# Compute the HEPData set indices assigned to one MC input
def dataset_indices_for_mc(assignment, dataset_sets):
    if assignment is None:
        return [index for index, dataset in enumerate(dataset_sets) if dataset.get("mc", True)]

    indices = [
        index
        for index, dataset in enumerate(dataset_sets)
        if assignment in steering.set_sample_names(dataset)
    ]
    if not indices:
        raise ValueError(f"iceplot: no dataset sets assigned to MC sample {assignment!r}")
    return indices


# Group dataset set indices by an optional shared plot directory
def dataset_plot_groups(dataset_sets):
    groups = []
    grouped_positions = {}
    for index, dataset in enumerate(dataset_sets):
        group = dataset.get("plot_group")
        output = group or dataset.get("plotname", dataset["name"])
        output = str(output).replace(" ", "_")
        output = f"{output}/" if output else ""
        if group is None:
            groups.append((output, [index]))
            continue
        if group not in grouped_positions:
            grouped_positions[group] = len(groups)
            groups.append((output, []))
        groups[grouped_positions[group]][1].append(index)
    return groups


# Apply one line style to every ordinary histogram record
def set_histogram_linestyle(histogram_sets, linestyle):
    if linestyle is None:
        return
    for histogram_set in histogram_sets:
        if histogram_set is None:
            continue
        for record in histogram_set.values():
            if "style" in record:
                record["style"] = {**record["style"], "ls": linestyle}


# Compute a plot copy with one shared channel color
def recolor_histogram_set(histogram_set, color):
    output = copy.deepcopy(histogram_set)
    for record in output.values():
        record["color"] = color
    return output


# Append one plot channel with palette colors and an optional fixed data color
def append_colored_plot_channel(
    grouped_mc, grouped_data, mc_sets, data_set, color_index, data_color=None, match_colors=False
):
    for histogram_set in mc_sets:
        grouped_mc.append(recolor_histogram_set(histogram_set, plots.colors(color_index)))
        if not match_colors:
            color_index += 1
    if data_set is not None:
        color = plots.colors(color_index) if match_colors else (data_color or plots.colors(color_index))
        grouped_data.append(recolor_histogram_set(data_set, color))
        if data_color is None and not match_colors:
            color_index += 1
    if match_colors and (mc_sets or data_set is not None):
        color_index += 1
    return color_index


# Compute the absolute HepMC3 path or glob pattern for one explicit source
def hepmc3_input_pattern(value: str, cdir: str) -> str:
    path = pathlib.Path(os.path.expanduser(value))
    if not path.is_absolute():
        path = pathlib.Path(cdir) / path
    return os.path.abspath(path)


# Resolve one logical HepMC3 source into a deterministic list of physical files
def resolve_hepmc3_source(value: str, cdir: str) -> list[str]:
    pattern = hepmc3_input_pattern(value=value, cdir=cdir)
    if os.path.isfile(pattern):
        return [pattern]

    files = sorted(path for path in glob.glob(pattern) if os.path.isfile(path))
    if not files:
        raise FileNotFoundError(f"HepMC3 input pattern '{pattern}' matched no files")
    return [os.path.abspath(path) for path in files]


# Resolve one analysis directory to its dataset steering card
def resolve_analysis_input(value: str, cdir: str) -> tuple[str, str]:
    analysis_dir = pathlib.Path(os.path.expanduser(value))
    if not analysis_dir.is_absolute():
        analysis_dir = pathlib.Path(cdir) / analysis_dir
    analysis_dir = pathlib.Path(os.path.abspath(analysis_dir))
    if not analysis_dir.is_dir():
        raise NotADirectoryError(f"Analysis directory '{analysis_dir}' not found")

    dataset_file = analysis_dir / "dataset.json"
    if not dataset_file.is_file():
        raise FileNotFoundError(f"Analysis dataset card '{dataset_file}' not found")
    return str(analysis_dir), str(dataset_file)


# Resolve a cut-definition file or built-in cuts module name
def resolve_cut_reference(
    value: str,
    cdir: str,
    dataset_path: str | None = None,
) -> str:
    return steering.resolve_python_reference(
        value,
        package="core.analysis.cuts",
        dataset_path=dataset_path,
        cdir=cdir,
    )


# Compute a filename-safe default tag for an input path or name
def input_tag(value: str) -> str:
    return pathlib.Path(value).stem


# Compute a compact tag for an input path, name or glob pattern
def hepmc3_source_tag(value: str) -> str:
    tag = input_tag(value)
    for token in ("*", "?", "[", "]"):
        tag = tag.replace(token, "")
    return tag.rstrip("._-") or "MC"


# Compute a bundle-aware tag for one icepack dataset card
def dataset_tag(value: str, cdir: str) -> str:
    path = pathlib.Path(value).resolve()
    if path.name != "dataset.json":
        return input_tag(value)
    try:
        relative = path.parent.relative_to(pathlib.Path(cdir).resolve() / "icepack")
        return "__".join(relative.parts)
    except ValueError:
        return path.parent.name


# Resolve one optional output path relative to the configured working directory
def resolve_output_path(value: str | None, cdir: str, default: str) -> str:
    path = pathlib.Path(default if value is None else os.path.expanduser(value))
    if not path.is_absolute():
        path = pathlib.Path(cdir) / path
    return str(path.resolve())


# Broadcast one cut definition or require one definition per HepMC3 input
def resolve_cli_cuts(cuts: list, hepmc3_count: int, base_dir: str) -> list:
    if hepmc3_count == 0:
        return []
    if len(cuts) == 1:
        cuts = cuts * hepmc3_count
    elif len(cuts) != hepmc3_count:
        raise ValueError(
            "iceplot: --cuts must provide one definition or one definition per HepMC3 input"
        )
    return [resolve_cut_reference(value=value, cdir=base_dir) for value in cuts]


# Parse command line arguments
def parse_args():
    """
    Parse input arguments
    """

    print(" ".join(sys.argv))

    parser = argparse.ArgumentParser(
        description=f"%(prog)s {__version__} ({__RELEASE__}) [{__AUTHOR__}]",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    cli.add_value(parser, "cdir", os.getcwd())
    cli.add_value(
        parser,
        "hepmc3",
        [],
        value_type=str,
        nargs="+",
        help="HepMC3 tags, file paths or quoted glob patterns",
    )
    cli.add_value(parser, "obs", "default", help="Standalone observable module")
    cli.add_value(
        parser,
        "cuts",
        ["default"],
        value_type=str,
        nargs="+",
        help="Cut modules or Python paths",
    )
    cli.add_value(parser, "pid", "[[211,-211]]", help="Final-state particle IDs")
    cli.add_value(parser, "mclabel", None, nargs="+", help="MC legend labels")
    cli.add_value(parser, "unit", None, help="Cross-section unit")
    cli.add_value(parser, "maxevents", None, value_type=int, help="Maximum events")
    cli.add_value(
        parser,
        "chunksize",
        "auto",
        value_type=parse_chunksize,
        help="Events per multiprocessing chunk",
    )
    cli.add_flag(parser, "chi2", help="Show chi2 in figure titles")
    cli.add_toggle(
        parser,
        "density",
        default=None,
        help="Normalize histograms to densities",
    )
    cli.add_flag(parser, "stack", help="Stack independent MC processes")
    cli.add_value(
        parser,
        "density-uncertainty",
        None,
        choices=["scaled", "shape"],
        help="Unit-density uncertainty propagation",
    )
    cli.add_value(
        parser,
        "ratio-uncertainty",
        None,
        choices=["combined", "separate", "numerator", "none"],
        help="Ratio-panel uncertainty convention",
    )
    cli.add_value(parser, "analysis", None, help="Analysis directory containing dataset.json")
    cli.add_value(
        parser,
        "mc-hepdata-sample",
        None,
        nargs="+",
        help="Dataset sample name assigned to each MC input",
    )
    cli.add_value(parser, "output", None, help="Output tag")
    cli.add_value(parser, "output-dir", None, help="Plot root directory")
    cli.add_value(parser, "report", None, help="Comparison report path")
    cli.add_flag(parser, "validate", help="Fail on unmet measurement criteria after saving outputs")
    cli.add_value(parser, "title", "", help="Figure title")
    cli.add_value(parser, "title_loc", "left", help="Figure title location")
    cli.add_value(
        parser,
        "legend",
        "inside",
        choices=["inside", "outside-right", "outside-top", "none"],
        help="Legend placement",
    )
    cli.add_value(parser, "legend-ncol", 1, help="Legend columns")
    cli.add_value(
        parser,
        "legend-fontsize",
        plots.DEFAULT_LEGEND_FONTSIZE,
        help="Legend font size",
    )
    cli.add_value(
        parser,
        "legend-bbox",
        None,
        value_type=float,
        nargs=2,
        help="Legend anchor in figure coordinates",
    )
    cli.add_value(
        parser,
        "datascale",
        None,
        value_type=float,
        nargs="+",
        help="Data scale factors",
    )
    cli.add_value(parser, "mcscale", None, nargs="+", help="MC scale factors")
    cli.add_value(parser, "xsmode", "auto", help="Cross-section normalization mode")
    cli.add_value(parser, "cores", 0, help="Maximum worker processes")
    cli.add_value(parser, "verbose", 1, help="Verbosity")
    args = parser.parse_args()
    args.unit_cli = args.unit is not None
    args.cores = max(1, multiprocessing.cpu_count() // 2) if (args.cores == 0) else args.cores
    if args.cores <= 0:
        raise ValueError("iceplot: --cores must be positive (use 0 for automatic)")

    args.pid = ast.literal_eval(args.pid)

    if args.xsmode not in ["auto", "header", "sample"]:
        raise ValueError("iceplot: Unknown xsmode. Use one of: auto, header, sample")

    if args.legend_ncol <= 0:
        raise ValueError("iceplot: --legend-ncol must be positive")
    if args.xsmode == "auto":
        cprint(
            "iceplot: Using 'xsmode' = 'auto'. Cross-section mode is resolved independently for each HepMC3 input.",
            "green",
        )
    elif args.xsmode == "sample":
        cprint(
            "iceplot: Using 'xsmode' = 'sample'. NOTE: this mode requires weighted (non-unit weight) events.",
            "green",
        )
    # Resolve path and module inputs

    args.cdir = os.path.abspath(os.path.expanduser(args.cdir))
    args.output_dir = resolve_output_path(
        value=args.output_dir,
        cdir=args.cdir,
        default="figs/iceplot",
    )
    args.report = (
        resolve_output_path(value=args.report, cdir=args.cdir, default=args.report)
        if args.report is not None
        else None
    )

    if (args.analysis is None) and (len(args.hepmc3) == 0):
        raise ValueError("iceplot: No HEPData or HepMC3 input given")
    args.hepmc3_sources = [
        resolve_hepmc3_source(value=value, cdir=args.cdir) for value in args.hepmc3
    ]
    args.hepmc3files = [
        files[0] if len(files) == 1 else hepmc3_input_pattern(value=value, cdir=args.cdir)
        for value, files in zip(args.hepmc3, args.hepmc3_sources, strict=True)
    ]
    args.analysis_dir = None
    args.dataset_file = None
    if args.analysis is not None:
        args.analysis_dir, args.dataset_file = resolve_analysis_input(
            value=args.analysis,
            cdir=args.cdir,
        )
    args.dataset_card = None
    if args.dataset_file is not None:
        args.dataset_card, args.dataset_file = steering.load_dataset(
            args.dataset_file,
            cdir=args.cdir,
        )
        plot = args.dataset_card["plot"]
        if args.unit is None:
            args.unit = plot.get("unit", "pb")
        args.stack = args.stack or plot["stack"]
        if args.density is None:
            args.density = plot["normalization"] == "unit_density"
        if args.density_uncertainty is None:
            args.density_uncertainty = (
                plot.get("density_uncertainty", "scaled") if args.density else "scaled"
            )
        if args.ratio_uncertainty is None:
            args.ratio_uncertainty = plot["ratio_uncertainty"]
    else:
        args.unit = args.unit or "pb"
        args.density = False if args.density is None else args.density
        args.density_uncertainty = args.density_uncertainty or "scaled"
        args.ratio_uncertainty = args.ratio_uncertainty or "combined"
    if args.stack and 0 < len(args.hepmc3) < 2:
        raise ValueError("iceplot: --stack requires at least two HepMC3 process sources")
    args.cut_references = (
        resolve_cli_cuts(cuts=args.cuts, hepmc3_count=len(args.hepmc3), base_dir=args.cdir)
        if args.analysis is None
        else []
    )
    # Create output folder name

    if args.output is None:
        if args.dataset_file is not None:
            tags = [dataset_tag(args.dataset_file, args.cdir)]
        else:
            tags = ["MC", *(hepmc3_source_tag(value) for value in args.hepmc3)]
        args.output = "__".join(tags)

    if args.density:
        args.output += "__[density]"
    if args.stack:
        args.output += "__[stack]"
    # Check input

    if args.mc_hepdata_sample is not None and args.analysis is None:
        raise ValueError("iceplot: --mc-hepdata-sample requires --analysis")

    if args.mcscale is not None and len(args.mcscale) != len(args.hepmc3):
        raise Exception("iceplot: Possible scale factors should be provided for each HepMC3 file")

    if args.mclabel is not None and len(args.mclabel) != len(args.hepmc3):
        raise Exception("iceplot: Possible text label should be provided for each HepMC3 file")

    if len(args.pid) != len(args.hepmc3) and (len(args.hepmc3) > 0 and (args.analysis is None)):
        args.pid = [args.pid[0]] * len(args.hepmc3)
        cprint(f"iceplot: Using same PID array for each HepMC3 file as {args.pid}", "red")

    return args


# Apply the selected cross-section unit to histograms and measured points
def change_scale(all_obs, args, units=None):
    """
    Change visualization scale
    """

    # These are w.r.t. to HepMC standard of picobarns
    unitset = plots.XS_UNITS

    selected = {key: units.get(key, args.unit) for key in all_obs} if units is not None else None
    requested = {args.unit} if selected is None else set(selected.values())
    unknown = requested.difference(unitset)
    if unknown:
        raise ValueError(
            f"Unknown unit '{sorted(unknown)[0]}'. Valid options: {', '.join(unitset.keys())}"
        )

    for key in all_obs:
        if plots.cross_section_observable(all_obs[key]):
            all_obs[key]["units"]["y"] = args.unit if selected is None else selected[key]

    if selected is None:
        return all_obs, unitset[args.unit]
    return all_obs, {key: unitset[unit] for key, unit in selected.items()}


# Combine one common scale with scalar or observable-specific unit scales
def apply_unit_scale(scale, unitscale):
    if isinstance(unitscale, dict):
        return {key: scale * value for key, value in unitscale.items()}
    return scale * unitscale


# Format a HepMC3 event weight for diagnostic text
def format_weight_diagnostic(value) -> str:
    """
    Format an event weight for warning diagnostics
    """
    if value is None:
        return "unavailable"
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return f"{float(value):.16E}"


# Resolve the effective cross-section mode from the requested mode and event weights
def resolve_xsmode(requested_xsmode: str, weight_summary: dict) -> tuple[str, str]:
    """
    Return the effective x-section mode and the reason for the choice
    """
    return readers.resolve_hepmc3_xsmode(
        requested_xsmode=requested_xsmode,
        weight_summary=weight_summary,
    )


# Combine per-file weight scans into one logical HepMC3 source summary
def combine_weight_summaries(summaries: list[dict]) -> dict:
    return readers.combine_hepmc3_weight_summaries(summaries)


# Scan one logical HepMC3 source while applying a global event limit
def scan_hepmc3_source(files: list[str], maxevents: int | None) -> tuple[dict, list[dict]]:
    remaining = maxevents
    summaries = []
    for filename in files:
        if remaining is not None and remaining <= 0:
            break
        summary = readers.read_hepmc3_weight_summary(
            filename=filename,
            maxevents=int(1e9) if remaining is None else remaining,
        )
        summaries.append(summary)
        if remaining is not None:
            remaining -= summary["nevents"]

    combined = combine_weight_summaries(summaries)
    combined["maxevents_reached"] = maxevents is not None and combined["nevents"] >= maxevents
    return combined, summaries


# Compute the standard weighted effective sample size divided by event count
def effective_sample_fraction(nevents: int, wsum: float, wsum2: float) -> float:
    if nevents <= 0 or wsum2 <= 0.0:
        return 0.0
    return min(wsum**2 / (nevents * wsum2), 1.0)


# Aggregate additive worker diagnostics into one logical HepMC3 source summary
def aggregate_source_diagnostics(chunks: list[dict], xsmode: str) -> dict:
    if not chunks:
        raise ValueError("iceplot: cannot aggregate an empty HepMC3 source")

    before = {
        key: sum(chunk["before"][key] for chunk in chunks)
        for key in ("nevents", "attempted_events", "wsum", "wsum2")
    }
    if xsmode == "sample":
        xsection_pb, xsection_pb_err = aggregate_sample_xsection(
            chunks=chunks,
            wsum=before["wsum"],
            wsum2=before["wsum2"],
        )
    else:
        xsection_pb, xsection_pb_err = aggregate_header_xsection(
            chunks=chunks,
            total_wsum=before["wsum"],
        )
    before.update(
        {
            "xsection_pb": xsection_pb,
            "xsection_pb_err": xsection_pb_err,
            "acceptance": 1.0,
            "event_acceptance": 1.0,
            "ess_fraction": effective_sample_fraction(
                before["nevents"], before["wsum"], before["wsum2"]
            ),
        }
    )

    after = []
    for set_index in range(len(chunks[0]["after"])):
        record = {
            key: sum(chunk["after"][set_index][key] for chunk in chunks)
            for key in ("nevents", "wsum", "wsum2")
        }
        acceptance = record["wsum"] / before["wsum"] if before["wsum"] != 0.0 else 0.0
        record.update(
            {
                "xsection_pb": xsection_pb * acceptance,
                "xsection_pb_err": aggregate_sample_xsection(chunks, record["wsum"], record["wsum2"])[1]
                if xsmode == "sample" else xsection_pb_err * abs(acceptance),
                "acceptance": acceptance,
                "event_acceptance": (
                    record["nevents"] / before["nevents"] if before["nevents"] > 0 else 0.0
                ),
                "ess_fraction": effective_sample_fraction(
                    record["nevents"], record["wsum"], record["wsum2"]
                ),
            }
        )
        after.append(record)
    return {"before": before, "after": after}


# Select the cut diagnostics routed to one MC input
def route_source_diagnostics(summary: dict, set_indices: list[int]) -> dict:
    routed = dict(summary)
    routed["after"] = [summary["after"][index] for index in set_indices]
    return routed


# Format one compact text table with explicit column delimiters
def format_text_table(headers: list[str], rows: list[tuple[str, ...]]) -> str:
    widths = [
        max([len(headers[index]), *(len(row[index]) for row in rows)])
        for index in range(len(headers))
    ]
    lines = [
        " | ".join(f"{header:<{width}}" for header, width in zip(headers, widths, strict=True)),
        "-+-".join("-" * width for width in widths),
    ]
    for row in rows:
        lines.append(
            " | ".join(
                f"{value:<{width}}" if index == 0 else f"{value:>{width}}"
                for index, (value, width) in enumerate(zip(row, widths, strict=True))
            )
        )
    return "\n".join(lines)


# Print one table and optionally save identical ASCII plus structured JSON outputs
def print_and_save_table(
    title: str,
    headers: list[str],
    rows: list[tuple[str, ...]],
    json_rows: list[dict],
    output_dir: str | pathlib.Path | None = None,
    table_name: str | None = None,
    metadata: dict | None = None,
) -> tuple[pathlib.Path, pathlib.Path] | None:
    if len(rows) != len(json_rows):
        raise ValueError("iceplot: table display and JSON row counts differ")

    table = format_text_table(headers=headers, rows=rows)
    file_title = title.lstrip("\n")
    cprint(title, "yellow")
    print(table)
    print()

    if output_dir is None and table_name is None:
        return None
    if output_dir is None or table_name is None:
        raise ValueError("iceplot: table output directory and name must be provided together")
    if pathlib.Path(table_name).name != table_name:
        raise ValueError("iceplot: table name must not contain path components")

    directory = pathlib.Path(output_dir)
    ensure_dir(directory)
    json_path = directory / f"{table_name}.json"
    text_path = directory / f"{table_name}.txt"
    payload = {
        "schema_version": 1,
        "title": file_title,
        "columns": headers,
        "metadata": {} if metadata is None else metadata,
        "rows": json_rows,
    }
    with json_path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        stream.write("\n")

    ascii_table = f"{file_title}\n{table}\n".encode("ascii", errors="backslashreplace").decode(
        "ascii"
    )
    text_path.write_text(ascii_table, encoding="ascii")
    cprint(f"Saving tables: {json_path}, {text_path}", "yellow")
    return json_path, text_path


# Format one source statistic with compact precision
def format_source_statistic(value: float, fraction: bool = False) -> str:
    return f"{float(value):.3f}" if fraction else f"{float(value):.3E}"


# Print one aligned before-cut and after-cut table for a logical HepMC3 source
def print_source_summary_table(
    source: str,
    summary: dict,
    names: list[str],
    output_dir: str | pathlib.Path | None = None,
    table_name: str | None = None,
) -> None:
    stages = [("Before cuts", summary["before"])]
    stages.extend(
        (
            f"After cuts [{name.rstrip('/')}]" if name.rstrip("/") else "After cuts",
            record,
        )
        for name, record in zip(names, summary["after"], strict=True)
    )
    headers = (
        "Stage",
        "events",
        "event acceptance",
        "mc_xs_pb",
        "mc_xs_pb_err",
        "cross section acceptance",
        "wsum",
        "wsum2",
        "ESS / N",
    )
    rows = [
        (
            stage,
            str(record["nevents"]),
            format_source_statistic(record["event_acceptance"], fraction=True),
            format_source_statistic(record["xsection_pb"]),
            format_source_statistic(record["xsection_pb_err"]),
            format_source_statistic(record["acceptance"], fraction=True),
            format_source_statistic(record["wsum"]),
            format_source_statistic(record["wsum2"]),
            format_source_statistic(record["ess_fraction"], fraction=True),
        )
        for stage, record in stages
    ]
    json_rows = [
        {
            "stage": stage,
            "events": int(record["nevents"]),
            "event_acceptance": float(record["event_acceptance"]),
            "mc_xs_pb": float(record["xsection_pb"]),
            "mc_xs_pb_err": float(record["xsection_pb_err"]),
            "cross_section_acceptance": float(record["acceptance"]),
            "wsum": float(record["wsum"]),
            "wsum2": float(record["wsum2"]),
            "ess_fraction": float(record["ess_fraction"]),
        }
        for stage, record in stages
    ]
    print_and_save_table(
        title=f"\niceplot: HepMC3 source summary [{source}]",
        headers=list(headers),
        rows=rows,
        json_rows=json_rows,
        output_dir=output_dir,
        table_name=table_name,
        metadata={"kind": "source_summary", "source": source},
    )


# Format one cross section or ratio with its uncertainty
def format_value_uncertainty(value: float | None, error: float | None) -> str:
    if value is None or error is None:
        return "undefined"
    return f"{value:.3E} +- {error:.3E}"


# Format one dimensionless ratio with fixed-point uncertainty
def format_ratio_uncertainty(value: float | None, error: float | None) -> str:
    if value is None or error is None:
        return "undefined"
    return f"{value:.3f} +- {error:.3f}"


# Compute the common valid-bin mask for one or more histograms
def histogram_valid_mask(*histograms):
    import numpy as np

    if not histograms:
        raise ValueError("iceplot: at least one histogram is required")
    size = len(np.asarray(histograms[0].counts_scaled, dtype=float))
    if len(histograms) > 1:
        reference = np.asarray(histograms[0].bins, dtype=float)
        if reference.ndim != 1 or reference.shape != (size + 1,):
            raise ValueError("iceplot: histogram bin edges have incompatible dimensions")
        if not np.all(np.isfinite(reference)) or np.any(np.diff(reference) <= 0.0):
            raise ValueError("iceplot: histogram bin edges must be finite and ordered")
        for histogram in histograms[1:]:
            bins = np.asarray(histogram.bins, dtype=float)
            if bins.shape != reference.shape or not np.array_equal(bins, reference):
                raise ValueError("iceplot: MC and data histograms use different bin edges")
    valid = np.ones(size, dtype=bool)
    for histogram in histograms:
        if len(np.asarray(histogram.counts_scaled, dtype=float)) != size:
            raise ValueError("iceplot: histograms have incompatible dimensions")
        local = np.asarray(
            getattr(histogram, "valid", np.ones(size, dtype=bool)),
            dtype=bool,
        )
        if local.shape != (size,):
            raise ValueError("iceplot: histogram valid-bin masks have incompatible dimensions")
        valid &= local
    return valid


# Compute the histogram integral over selected bins
def histogram_integral(histogram, valid=None):
    import numpy as np

    values = np.asarray(histogram.counts_scaled, dtype=float)
    weights = np.asarray(histogram.measure, dtype=float)
    mask = histogram_valid_mask(histogram) if valid is None else np.asarray(valid, dtype=bool)
    if values.shape != weights.shape or mask.shape != values.shape:
        raise ValueError("iceplot: histogram integral arrays have incompatible dimensions")
    return float(np.sum(np.where(mask, values * weights, 0.0)))


# Compute the integral uncertainty with source correlations
def histogram_integral_error(histogram, uncertainties=None, valid=None):
    import numpy as np

    covariance = histogram_covariance(histogram, uncertainties)
    weights = np.asarray(histogram.measure, dtype=float)
    mask = histogram_valid_mask(histogram) if valid is None else np.asarray(valid, dtype=bool)
    if mask.shape != weights.shape:
        raise ValueError("iceplot: histogram integral mask has incompatible dimensions")
    weights = np.where(mask, weights, 0.0)
    return float(uncertainty.covariance_projection_error(covariance, weights))


# Compute a covariance-aware MC and data chi-square when sources are available
def histogram_chi2(mc_histogram, data_histogram, data_uncertainties=None):
    import numpy as np

    mc_values = np.asarray(mc_histogram.counts_scaled, dtype=float)
    data_values = np.asarray(data_histogram.counts_scaled, dtype=float)
    mc_errors = np.asarray(mc_histogram.errs_scaled, dtype=float)
    if mc_values.shape != data_values.shape or mc_values.shape != mc_errors.shape:
        raise ValueError("iceplot: MC and data histograms have incompatible dimensions")
    valid = np.asarray(
        getattr(mc_histogram, "valid", np.ones(len(mc_values), dtype=bool)),
        dtype=bool,
    ) & np.asarray(
        getattr(data_histogram, "valid", np.ones(len(data_values), dtype=bool)),
        dtype=bool,
    )
    result = objective.correlated_chi2_arrays(
        mc_prediction=mc_values,
        data_values=data_values,
        mc_stat_uncertainty=mc_errors,
        fit_weights=valid.astype(float),
        data_total_covariance=histogram_covariance(data_histogram, data_uncertainties),
        mc_covariance=mc_histogram.covariance_scaled,
    )
    if not result["valid"]:
        return 0.0, 0
    return float(result["chi2"]), int(result["rank"])


# Print histogram-integral summaries using the selected normalization
def print_mc_data_comparison_tables(
    source: str | None,
    mc_sets: list[dict] | None,
    data_sets: list[dict] | None,
    names: list[str],
    unit: str,
    normalization: str = "cross_section",
    output_dir: str | pathlib.Path | None = None,
    table_prefix: str | None = None,
) -> None:
    if normalization not in {"cross_section", "unit_density"}:
        raise ValueError("iceplot: unknown histogram summary normalization")
    if mc_sets is None and data_sets is None:
        return
    density = normalization == "unit_density"
    mc_key = "mc_integral" if density else "mc_xs"
    data_key = "data_integral" if density else "data_xs"
    set_count = len(mc_sets) if mc_sets is not None else len(data_sets)
    if mc_sets is not None and data_sets is not None and len(data_sets) != set_count:
        raise ValueError("iceplot: MC and data histogram set counts differ")

    for set_index in range(set_count):
        mc_set = None if mc_sets is None else mc_sets[set_index]
        data_set = None if data_sets is None else data_sets[set_index]
        has_mc = mc_set is not None
        has_data = data_set is not None
        observables = list(
            dict.fromkeys(
                [
                    *(mc_set.keys() if mc_set is not None else ()),
                    *(data_set.keys() if data_set is not None else ()),
                ]
            )
        )
        observable_units = {}
        for observable in observables:
            mc_record = None if mc_set is None else mc_set.get(observable)
            data_record = None if data_set is None else data_set.get(observable)
            reference = mc_record or data_record or {}
            observable_units[observable] = (
                "1"
                if density
                else reference.get("obs", {}).get("units", {}).get("y", unit)
            )
        mixed_units = len(set(observable_units.values())) > 1
        rows = []
        json_rows = []
        for observable in observables:
            mc_record = None if mc_set is None else mc_set.get(observable)
            data_record = None if data_set is None else data_set.get(observable)
            reference_record = mc_record or data_record
            if (
                reference_record is not None
                and reference_record.get("obs", {}).get("kind", "histogram") == "point"
            ):
                continue
            mc_histogram = mc_record.get("hdata") if isinstance(mc_record, dict) else None
            data_histogram = data_record.get("hdata") if isinstance(data_record, dict) else None
            if mc_histogram is None and data_histogram is None:
                continue
            common_valid = (
                histogram_valid_mask(mc_histogram, data_histogram)
                if mc_histogram is not None and data_histogram is not None
                else None
            )

            row = [observable]
            if mixed_units:
                row.append(observable_units[observable])
            mc_values = None
            data_values = None
            if mc_histogram is not None:
                mc_values = (
                    histogram_integral(mc_histogram, valid=common_valid),
                    histogram_integral_error(mc_histogram, valid=common_valid),
                )
                row.append(format_value_uncertainty(*mc_values))
            elif has_mc:
                row.append("undefined")
            if data_histogram is not None:
                data_values = (
                    histogram_integral(data_histogram, valid=common_valid),
                    float(
                        histogram_integral_error(
                            data_histogram,
                            data_record.get("uncertainties", []),
                            valid=common_valid,
                        )
                    ),
                )
                row.append(format_value_uncertainty(*data_values))
            elif has_data:
                row.append("undefined")
            if has_mc and has_data:
                ratio = (
                    (None, None)
                    if mc_values is None or data_values is None
                    else cross_section_ratio(
                        numerator=mc_values[0],
                        numerator_error=mc_values[1],
                        denominator=data_values[0],
                        denominator_error=data_values[1],
                    )
                )
                row.append(format_ratio_uncertainty(*ratio))
            rows.append(tuple(row))
            json_row = {"observable": observable, "unit": observable_units[observable]}
            if has_mc:
                json_row[mc_key] = (
                    None if mc_values is None else {"value": mc_values[0], "error": mc_values[1]}
                )
            if has_data:
                json_row[data_key] = (
                    None
                    if data_values is None
                    else {"value": data_values[0], "error": data_values[1]}
                )
            if has_mc and has_data:
                json_row["mc_to_data"] = (
                    None if ratio[0] is None else {"value": ratio[0], "error": ratio[1]}
                )
            json_rows.append(json_row)

        if not rows:
            continue
        set_name = names[set_index].rstrip("/")
        suffix = f", set={set_name}" if set_name else ""
        headers = ["Observable"]
        if mixed_units:
            headers.append("Unit")
        display_unit = unit if mixed_units else next(iter(observable_units.values()), unit)
        if has_mc:
            headers.append("mc_integral +-" if density else f"mc_xs_{display_unit} +-")
        if has_data:
            headers.append("data_integral +-" if density else f"data_xs_{display_unit} +-")
        if has_mc and has_data:
            headers.append("mc_integral / data_integral +-" if density else "mc_xs / data_xs +-")
        data_label = set_name or "HEPData"
        source_label = f" [{source}{suffix}]" if source is not None else f" [{data_label}]"
        table_name = None if table_prefix is None else f"{table_prefix}_set_{set_index:03d}"
        summary_name = "unit-density" if density else "cross-section"
        print_and_save_table(
            title=f"iceplot: Histogram {summary_name} summary{source_label}",
            headers=headers,
            rows=rows,
            json_rows=json_rows,
            output_dir=output_dir,
            table_name=table_name,
            metadata={
                "kind": f"histogram_{normalization}_summary",
                "normalization": normalization,
                "source": source,
                "set": set_name,
                "unit": "1" if density else ("mixed" if mixed_units else display_unit),
            },
        )


# Build summed MC process sets for stack diagnostics
def build_summed_mc_process_sets(
    args,
    mclist: list,
    names: list[str],
    dataset_card: dict | None,
) -> list:
    stacked_sets = []
    for set_index in range(len(names)):
        indices = list(
            mc_indices_for_dataset(
                assignments=args.mc_hepdata_sample,
                dataset=dataset_card["sets"][set_index] if dataset_card is not None else {},
                input_count=len(mclist),
            )
        )
        if not indices:
            stacked_sets.append(None)
            continue

        reference = mclist[indices[0]][set_index]
        if any(set(mclist[index][set_index]) != set(reference) for index in indices):
            raise ValueError(
                f"iceplot: process stack inputs disagree on observables in set {set_index}"
            )
        stacked = {}
        for observable in reference:
            records = [mclist[index][set_index][observable] for index in indices]
            _, total = plots.build_filled_process_stack_records(
                records,
                density=args.density,
                density_uncertainty=args.density_uncertainty,
            )
            stacked[observable] = total
        stacked_sets.append(stacked)
    return stacked_sets


# Compute a concise source label for the selected x-section mode
def format_xsmode_source(
    requested_xsmode: str,
    resolved_xsmode: str,
    reason: str | None = None,
) -> str:
    """
    Format the source of the effective x-section mode selection
    """
    if requested_xsmode == "auto":
        if reason is not None and "LHE event weights" in reason:
            return "auto: LHE event weights"
        if reason is not None and "maximum weight overflow" in reason:
            return "auto: maximum weight overflow ratios"
        if resolved_xsmode == "header":
            return "auto: constant weights"
        if resolved_xsmode == "sample":
            return "auto: variable weights"
        return "auto"
    return "explicit"


# Print the effective cross-section mode selected for one HepMC3 input
def print_xsmode_resolution(
    hepmc3file: str,
    requested_xsmode: str,
    resolved_xsmode: str,
    reason: str | None = None,
) -> None:
    """
    Print the effective cross-section normalization mode used for one HepMC3 input
    """
    filename = os.path.basename(hepmc3file)
    source = format_xsmode_source(
        requested_xsmode=requested_xsmode,
        resolved_xsmode=resolved_xsmode,
        reason=reason,
    )
    cprint(f"iceplot: MC cross-section mode [{filename}]: {resolved_xsmode} ({source})", "green")


# Compute matplotlib legend keyword arguments from parsed CLI arguments
def legend_properties_from_args(args) -> dict:
    """
    Return legend keyword arguments for plot.superplot
    """
    return {
        "fontsize": args.legend_fontsize,
        "ncol": args.legend_ncol,
    }


# Print a prominent warning for header x-section mode with variable event weights
def warn_header_xsmode_weighted_events(hepmc3file: str, weight_summary: dict, xsmode: str) -> None:
    """
    Warn when header cross-section normalization is used on non-constant event weights
    """
    if (
        xsmode != "header"
        or not weight_summary.get("is_nonconstant", False)
        or weight_summary.get("maximum_weight_overflow_count", 0) > 0
    ):
        return

    first_event = weight_summary.get("first_nonconstant_event")
    if first_event is None:
        first_event_text = "not isolated"
    else:
        first_event_text = (
            f"{first_event} "
            f"(reference={format_weight_diagnostic(weight_summary.get('reference_weight'))}, "
            f"event_weight={format_weight_diagnostic(weight_summary.get('first_nonconstant_weight'))})"
        )

    lines = [
        "",
        "!" * 96,
        "!!! iceplot WARNING: --xsmode header with non-constant HepMC3 event weights !!!",
        f"File: {hepmc3file}",
        f"Events scanned: {weight_summary.get('nevents', 0)}",
        (
            f"Finite first-weight range: "
            f"[{format_weight_diagnostic(weight_summary.get('min_weight'))}, "
            f"{format_weight_diagnostic(weight_summary.get('max_weight'))}]"
        ),
        f"First non-constant event: {first_event_text}",
        (
            f"Missing weights: {weight_summary.get('missing_weight_count', 0)} | "
            f"non-finite weights: {weight_summary.get('nonfinite_weight_count', 0)}"
        ),
        "Header mode uses the HepMC3 cross-section header for normalization while the sample has variable event weights.",
        "Use --xsmode auto or --xsmode sample for weighted events, or regenerate/use a sample with constant event weights.",
        "!" * 96,
        "",
    ]

    for line in lines:
        cprint(line, "red", attrs=["bold"])


# Print a prominent warning for sample mode with maximum weight overflow ratios
def warn_sample_xsmode_overflow_events(hepmc3file: str, weight_summary: dict, xsmode: str) -> None:
    """
    Warn when sample cross-section normalization is used on dimensionless overflow ratios
    """
    overflow_count = weight_summary.get("maximum_weight_overflow_count", 0)
    if xsmode != "sample" or overflow_count == 0:
        return

    lines = [
        "",
        "!" * 96,
        "!!! iceplot WARNING: --xsmode sample with maximum weight overflow ratios !!!",
        f"File: {hepmc3file}",
        f"Events scanned: {weight_summary.get('nevents', 0)}",
        f"Maximum weight overflow events: {overflow_count}",
        "These event weights are dimensionless ratios and do not carry the sample cross section.",
        "Use --xsmode auto or --xsmode header so the HepMC3 cross-section header sets the normalization.",
        "!" * 96,
        "",
    ]

    for line in lines:
        cprint(line, "red", attrs=["bold"])


# Print a prominent warning for sample x-section mode with constant event weights
def warn_sample_xsmode_constant_weights(hepmc3file: str, weight_summary: dict, xsmode: str) -> None:
    """
    Warn when event-weight cross-section normalization is used on constant event weights
    """
    if xsmode != "sample" or not weight_summary.get("is_constant", False):
        return

    lines = [
        "",
        "!" * 96,
        "!!! iceplot WARNING: --xsmode sample with constant HepMC3 event weights !!!",
        f"File: {hepmc3file}",
        f"Events scanned: {weight_summary.get('nevents', 0)}",
        f"Constant first weight: {format_weight_diagnostic(weight_summary.get('reference_weight'))}",
        f"Max absolute deviation: {format_weight_diagnostic(weight_summary.get('max_abs_deviation'))}",
        "Sample mode computes the MC cross section from event weights, but the scanned event weights are constant.",
        "This usually indicates an unweighted/header-normalized sample; --xsmode auto or --xsmode header is likely the intended mode.",
        "Check that the event weights really carry the cross-section normalization before trusting --xsmode sample.",
        "!" * 96,
        "",
    ]

    for line in lines:
        cprint(line, "red", attrs=["bold"])


# Compose the CLI, set-selection and optional chi-square title lines
def compose_plot_title(cli_title, set_title=None, chi2=None, ndf=None):
    parts = [value for value in (cli_title, set_title) if value]
    if chi2 is not None:
        if ndf is None:
            raise ValueError("iceplot: chi-square plot title requires ndf")
        reduced = f"{chi2 / ndf:0.1f}" if ndf > 0 else "undefined"
        parts.append(f"$\\chi^2$ / ndf = {chi2:0.1f} / {ndf:0.0f} = {reduced}")
    return "\n".join(parts)


# Combine the sample-wide visualization scale with each dataset-set scale
def effective_mc_scales(base_scale, set_scales):
    return [base_scale * set_scale for set_scale in set_scales]


# Select the configured MC observables without loading a data reference
def select_mc_only_observables(dataset_set, all_obs):
    return steering.select_histogram_observables(dataset_set, all_obs)


# Check whether one independent set or shared plot group uses a ratio panel
def plot_group_ratio_enabled(dataset_card, set_indices):
    default = len(set_indices) == 1
    if dataset_card is None:
        return default
    sets = dataset_card.get("sets")
    if sets is not None and all(sets[index].get("mc", True) is False for index in set_indices):
        return False
    return dataset_card.get("plot", {}).get("ratio_plot", default)


# Create and save one configured ROC figure
def create_roc_figure(mc_obj, OBS, name, title, args):
    if len(mc_obj) < 2:
        raise ValueError("iceplot: ROC plotting requires at least two MC samples")

    fig, _ = plots.rocplot(
        mc_obj,
        title=title,
        title_loc=args.title_loc,
        legend_position=args.legend,
        legend_bbox=args.legend_bbox,
        legend_properties=legend_properties_from_args(args),
    )
    fullpath = os.path.join(args.output_dir, args.output, readers.clean_filename(name), "roc")
    ensure_dir(pathlib.Path(fullpath))
    filename = os.path.join(fullpath, f"rocplot__{OBS}.pdf")
    fig.savefig(filename, bbox_inches="tight")
    cprint(f"Saving figure: {filename}\n", "yellow")
    plt.close(fig)


# Create and save one histogram or ROC observable figure
def create_figures_parallel(param):
    """
    Parallel processing wrapper
    """

    # Take in parameters
    data = param["data"]
    mc = param["mc"]
    OBS = param["OBS"]
    name = param["name"]
    set_title = param["title"]
    args = param["args"]
    ratio_plot = param.get("ratio_plot", True)

    # Create
    data_obj = [data[i][OBS] for i in range(len(data))]
    mc_obj = [mc[i][OBS] for i in range(len(mc))]
    hist_objs = mc_obj + data_obj

    observable = hist_objs[0]["obs"]
    if observable.get("kind", "histogram") == "roc":
        if args.stack:
            raise ValueError("iceplot: --stack does not support ROC observables")
        if len(data_obj) > 0:
            raise ValueError("iceplot: ROC observables currently compare MC samples only")
        create_roc_figure(
            mc_obj=mc_obj,
            OBS=OBS,
            name=name,
            title=compose_plot_title(args.title, set_title),
            args=args,
        )
        return

    stack_total = None
    ratio_records = None
    legend_order = None
    stack_active = args.stack and bool(mc_obj)
    if stack_active:
        stacked_mc, stack_total = plots.build_filled_process_stack_records(
            mc_obj,
            density=args.density,
            density_uncertainty=args.density_uncertainty,
        )
        data_obj = plots.foreground_histogram_records(data_obj, stacked_mc + [stack_total])
        hist_objs = stacked_mc + data_obj
        ratio_records = data_obj + [stack_total] if data_obj else None
        legend_order = [record["label"] for record in data_obj + mc_obj]
    elif param.get("mc_first", False):
        mc_obj = plots.foreground_histogram_records(mc_obj, data_obj)
        hist_objs = mc_obj + data_obj
        ratio_records = data_obj + mc_obj if data_obj else None
        legend_order = [record["label"] for record in mc_obj + data_obj]
    else:
        data_obj = plots.foreground_histogram_records(data_obj, mc_obj)
        hist_objs = mc_obj + data_obj
        if data_obj:
            ratio_records = data_obj + mc_obj
            legend_order = [record["label"] for record in data_obj + mc_obj]

    for yscale in ["linear", "log"]:
        fig, ax = plots.superplot(
            hist_objs,
            ratio_plot=ratio_plot and (not stack_active or bool(data_obj)),
            ratio_data=ratio_records,
            uncertainty_record=stack_total,
            ratio_uncertainty=args.ratio_uncertainty,
            legend_order=legend_order,
            yscale=yscale,
            verbose=False,
            legend_position=args.legend,
            legend_bbox=args.legend_bbox,
            legend_properties=legend_properties_from_args(args),
        )

        ### Compute chi2 for this observable
        # ndf is the retained covariance rank
        if len(data_obj) == 1 and len(mc_obj) == 1 and args.chi2:
            i, j = 0, 0  # Use the first pair
            mc_histogram = stack_total["hdata"] if stack_total is not None else mc_obj[i]["hdata"]
            chi2, nbins = histogram_chi2(
                mc_histogram=mc_histogram,
                data_histogram=data_obj[j]["hdata"],
                data_uncertainties=data_obj[j].get("uncertainties", []),
            )

            title = compose_plot_title(args.title, set_title, chi2=chi2, ndf=nbins)
            ax[0].set_title(title, loc=args.title_loc, fontsize=6)
        else:
            title = compose_plot_title(args.title, set_title)
            ax[0].set_title(title, loc=args.title_loc, fontsize=6)

        # Create path and save
        fullpath = os.path.join(args.output_dir, args.output, readers.clean_filename(name), yscale)
        ensure_dir(pathlib.Path(fullpath))
        filename = os.path.join(fullpath, f"hplot__{OBS}.pdf")

        fig.savefig(filename, bbox_inches="tight")
        cprint(f"Saving figure: {filename}\n", "yellow")
        plt.close()


# Compute finite histogram values used by the validation report
def histogram_summary(histogram, allow_empty=False, uncertainties=None, valid=None):
    import numpy as np

    integral = histogram_integral(histogram, valid=valid)
    error = float(histogram_integral_error(histogram, uncertainties, valid=valid))
    counts = np.asarray(histogram.counts_scaled, dtype=float)
    errors = np.asarray(histogram.errs_scaled, dtype=float)
    mask = histogram_valid_mask(histogram) if valid is None else np.asarray(valid, dtype=bool)
    if counts.shape != errors.shape or mask.shape != counts.shape:
        raise ValueError("iceplot: histogram summary arrays have incompatible dimensions")
    empty = not np.any(mask & ((counts != 0.0) | (errors != 0.0)))
    if not math.isfinite(integral) or not math.isfinite(error):
        raise ValueError("iceplot: non-finite histogram integral in validation report")
    if not np.isfinite(counts).all() or not np.isfinite(errors).all():
        raise ValueError("iceplot: non-finite histogram bins in validation report")
    if error < 0.0 or np.any(errors < 0.0):
        raise ValueError("iceplot: negative histogram uncertainty in validation report")
    if integral < 0.0 or (empty and not allow_empty):
        raise ValueError("iceplot: invalid histogram integral in validation report")
    summary = {
        "integral": integral,
        "integral_error": error,
        "empty": empty,
    }
    summary.update(histogram_mc_precision(counts, errors, mask))
    summary.update(uncertainty.mc_statistics(histogram, valid=mask))
    return summary


# Compute finite ROC sample values used by the validation report
def roc_summary(rocdata):
    import numpy as np

    scores = np.asarray(rocdata["scores"], dtype=float)
    weights = np.asarray(rocdata["weights"], dtype=float)
    if scores.ndim != 2 or scores.shape[0] == 0:
        raise ValueError("iceplot: empty or malformed ROC scores in validation report")
    if weights.shape != (scores.shape[0],):
        raise ValueError("iceplot: ROC score and weight counts differ in validation report")
    if not np.isfinite(scores).all() or not np.isfinite(weights).all():
        raise ValueError("iceplot: non-finite ROC values in validation report")
    return {
        "events": int(scores.shape[0]),
        "weight_sum": float(np.sum(weights)),
    }


# Compute literature-comparison metrics for one MC and data histogram pair
def comparison_metrics(mc_histogram, data_histogram, data_uncertainties=None):
    import numpy as np

    common_valid = histogram_valid_mask(mc_histogram, data_histogram)
    mc_summary = histogram_summary(mc_histogram, allow_empty=True, valid=common_valid)
    data_summary = histogram_summary(
        data_histogram,
        uncertainties=data_uncertainties,
        valid=common_valid,
    )
    chi2, ndf = histogram_chi2(
        mc_histogram=mc_histogram,
        data_histogram=data_histogram,
        data_uncertainties=data_uncertainties,
    )
    chi2_available = ndf > 0 and math.isfinite(float(chi2))
    combined_error = math.hypot(
        mc_summary["integral_error"],
        data_summary["integral_error"],
    )
    integral_pull = (
        (mc_summary["integral"] - data_summary["integral"]) / combined_error
        if combined_error > 0.0
        else None
    )
    shape_l1 = None
    if mc_summary["integral"] > 0.0 and data_summary["integral"] > 0.0:
        binwidth = np.asarray(data_histogram.measure, dtype=float)
        mc_density = np.where(
            common_valid,
            np.asarray(mc_histogram.counts_scaled, dtype=float) / mc_summary["integral"],
            0.0,
        )
        data_density = np.where(
            common_valid,
            np.asarray(data_histogram.counts_scaled, dtype=float) / data_summary["integral"],
            0.0,
        )
        shape_l1 = float(np.sum(np.abs(mc_density - data_density) * binwidth))
        if not math.isfinite(shape_l1):
            raise ValueError("iceplot: non-finite shape distance in validation report")
    return {
        "integral": mc_summary["integral"],
        "integral_error": mc_summary["integral_error"],
        "empty": mc_summary["empty"],
        "mc_max_rel_uncertainty": mc_summary["mc_max_rel_uncertainty"],
        "mc_zero_bin_count": mc_summary["mc_zero_bin_count"],
        "mc_empty_bin_count": mc_summary["mc_empty_bin_count"],
        "mc_empty_bins": mc_summary["mc_empty_bins"],
        "mc_min_neff": mc_summary["mc_min_neff"],
        "mc_valid_bin_count": mc_summary["mc_valid_bin_count"],
        "data_valid_bin_count": int(histogram_valid_mask(data_histogram).sum()),
        "mc_data_error_ratio": uncertainty.covariance_error_ratio(
            histogram_covariance(mc_histogram)[np.ix_(common_valid, common_valid)],
            histogram_covariance(data_histogram, data_uncertainties)[np.ix_(common_valid, common_valid)],
        ),
        "comparison_status": "empty_mc" if mc_summary["empty"] else "compared",
        "data_integral": data_summary["integral"],
        "data_integral_error": data_summary["integral_error"],
        "integral_pull": integral_pull,
        "chi2": float(chi2) if chi2_available else None,
        "ndf": int(ndf),
        "chi2_ndf": float(chi2 / ndf) if chi2_available else None,
        "chi2_status": "available" if chi2_available else "unavailable",
        "shape_l1": shape_l1,
    }


# Build the summed prediction row used by a stacked validation report
def stacked_validation_sample(args, mclist, mc_indices, set_index, observable, data_record):
    records = [mclist[index][set_index][observable] for index in mc_indices]
    _, total = plots.build_filled_process_stack_records(
        records,
        density=args.density,
        density_uncertainty=args.density_uncertainty,
    )
    row = {
        "label": total["label"],
        "inputs": [args.hepmc3files[index] for index in mc_indices],
        "input_files": [args.hepmc3_sources[index] for index in mc_indices],
        "components": [
            {
                "input": args.hepmc3files[index],
                "input_files": args.hepmc3_sources[index],
                "label": mclist[index][set_index][observable]["label"],
            }
            for index in mc_indices
        ],
    }
    if data_record is None:
        row.update(histogram_summary(total["hdata"], allow_empty=True))
    else:
        row.update(
            comparison_metrics(
                total["hdata"],
                data_record["hdata"],
                data_record.get("uncertainties", []),
            )
        )
    return row


# Compute the input index and label of an optional MC reference sample
def mc_validation_reference(dataset_card, mc_indices):
    if dataset_card is None:
        return None
    reference_name = dataset_card.get("validation", {}).get("mc_reference")
    if reference_name is None:
        return None
    sample_index = next(
        index
        for index, sample in enumerate(dataset_card["samples"])
        if sample["name"] == reference_name
    )
    if sample_index not in mc_indices:
        raise ValueError(
            f"iceplot: MC reference sample '{reference_name}' is not assigned to this set"
        )
    return sample_index, dataset_card["samples"][sample_index]["label"]


# Compute post-cut source statistics for one report sample and dataset set
def source_validation_metrics(diagnostics, set_index):
    if diagnostics is None:
        return {}
    after = diagnostics.get("after") if isinstance(diagnostics, dict) else None
    if not isinstance(after, list) or not 0 <= set_index < len(after):
        raise ValueError("iceplot: source diagnostics do not match validation sets")
    record = after[set_index]
    events = int(record["nevents"])
    ess_fraction = float(record["ess_fraction"])
    if events < 0 or not math.isfinite(ess_fraction) or not 0.0 <= ess_fraction <= 1.0:
        raise ValueError("iceplot: invalid source effective sample diagnostics")
    return {
        "mc_selected_events": events,
        "mc_ess_fraction": ess_fraction,
        "mc_effective_events": ess_fraction * events,
    }


# Build one machine-readable report from all iceplot histogram objects
def build_validation_report(
    args,
    mclist,
    data,
    names,
    dataset_card,
    source_diagnostics=None,
):
    if source_diagnostics is None:
        source_diagnostics = [None] * len(mclist)
    if len(source_diagnostics) != len(mclist):
        raise ValueError("iceplot: source diagnostics do not match MC inputs")
    sets = []
    nsets = len(names)
    for set_index in range(nsets):
        dataset_set = dataset_card["sets"][set_index] if dataset_card is not None else {}
        assignments = (
            mc_indices_for_dataset(
                assignments=args.mc_hepdata_sample,
                dataset=dataset_set,
                input_count=len(mclist),
            )
            if dataset_card is not None
            else range(len(mclist))
        )
        mc_indices = list(assignments)
        mc_reference = mc_validation_reference(dataset_card, mc_indices)
        observables = (
            list(mclist[mc_indices[0]][set_index]) if mc_indices else list(data[set_index] or {})
        )
        observable_rows = []
        for observable in observables:
            first_object = (
                mclist[mc_indices[0]][set_index][observable]
                if mc_indices
                else data[set_index][observable]
            )
            if "rocdata" in first_object:
                if args.stack:
                    raise ValueError("iceplot: process stacking does not support ROC reports")
                observable_rows.append(
                    {
                        "observable": observable,
                        "kind": "roc",
                        "samples": [
                            {
                                "input": args.hepmc3files[mc_index],
                                "input_files": args.hepmc3_sources[mc_index],
                                "label": mclist[mc_index][set_index][observable]["label"],
                                **roc_summary(mclist[mc_index][set_index][observable]["rocdata"]),
                            }
                            for mc_index in mc_indices
                        ],
                    }
                )
                continue
            data_record = (
                data[set_index][observable]
                if data is not None and data[set_index] is not None
                else None
            )
            if not mc_indices:
                samples = []
            elif args.stack:
                samples = [
                    stacked_validation_sample(
                        args=args,
                        mclist=mclist,
                        mc_indices=mc_indices,
                        set_index=set_index,
                        observable=observable,
                        data_record=data_record,
                    )
                ]
            else:
                samples = []
                for mc_index in mc_indices:
                    mc_object = mclist[mc_index][set_index][observable]
                    row = {
                        "input": args.hepmc3files[mc_index],
                        "input_files": args.hepmc3_sources[mc_index],
                        "label": mc_object["label"],
                    }
                    if data_record is not None:
                        row.update(
                            comparison_metrics(
                                mc_object["hdata"],
                                data_record["hdata"],
                                data_record.get("uncertainties", []),
                            )
                        )
                    elif mc_reference is not None and mc_index != mc_reference[0]:
                        reference_object = mclist[mc_reference[0]][set_index][observable]
                        row.update(
                            comparison_metrics(
                                mc_object["hdata"],
                                reference_object["hdata"],
                            )
                        )
                        row["comparison_kind"] = "mc_reference"
                        row["comparison_reference"] = mc_reference[1]
                    elif mc_reference is not None:
                        row.update(histogram_summary(mc_object["hdata"], allow_empty=True))
                        row["comparison_status"] = "reference"
                        row["comparison_kind"] = "mc_reference"
                        row["comparison_reference"] = mc_reference[1]
                    else:
                        row.update(histogram_summary(mc_object["hdata"], allow_empty=True))
                    row.update(
                        source_validation_metrics(
                            source_diagnostics[mc_index],
                            set_index,
                        )
                    )
                    samples.append(row)
            if not dataset_set.get("data", True) and mc_reference is None:
                for sample in samples:
                    sample["comparison_status"] = "prediction"
            observable_rows.append(
                {
                    "observable": observable,
                    "kind": first_object.get("obs", {}).get("kind", "histogram"),
                    "samples": samples,
                }
            )
        plot_name, _ = dataset_plot_groups([{"name": names[set_index], **dataset_set}])[0]
        sets.append(
            {
                "name": readers.clean_filename(names[set_index].rstrip("/")),
                "plotname": readers.clean_filename(plot_name.rstrip("/")),
                "title": dataset_set.get("title"),
                "region": dataset_set.get("region"),
                "mc": dataset_set.get("mc", True),
                "mc_scale": float(dataset_set.get("mc_scale", 1.0)),
                "observables": observable_rows,
            }
        )
    return {
        "schema_version": 1,
        "output": args.output,
        "plot_directory": str(pathlib.Path(args.output_dir) / args.output),
        "hepdata": args.dataset_file,
        "normalization": "unit_density" if args.density else "cross_section",
        "reference_normalization": dataset_card["plot"]["normalization"] if dataset_card and "plot" in dataset_card
        else "unit_density" if args.density else "cross_section",
        "density_uncertainty": args.density_uncertainty if args.density else None,
        "ratio_uncertainty": args.ratio_uncertainty,
        "stack": bool(args.stack),
        "validation": copy.deepcopy(dataset_card.get("validation", {}))
        if dataset_card is not None
        else {},
        "sets": sets,
    }


# Format one optional floating-point table field
def table_value(value):
    return "" if value is None else f"{float(value):.8g}"


# Write JSON and tab-separated iceplot validation tables
def write_validation_report(report, path):
    report["measurement"] = measurement_assessment(report)
    report_path = pathlib.Path(path)
    ensure_dir(report_path.parent)
    with report_path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")

    columns = (
        "set",
        "observable",
        "sample",
        "comparison_status",
        "mc_empty",
        "mc_integral",
        "mc_error",
        "mc_max_rel_uncertainty",
        "mc_zero_bin_count",
        "mc_effective_events",
        "mc_ess_fraction",
        "data_integral",
        "data_error",
        "integral_pull",
        "chi2",
        "ndf",
        "chi2_ndf",
        "shape_l1",
        "mc_data_error_ratio",
        "measurement_status",
        "p_chi2",
        "p_integral",
    )
    lines = ["\t".join(columns)]
    for dataset in report["sets"]:
        for observable in dataset["observables"]:
            for sample in observable["samples"]:
                values = (
                    dataset["name"],
                    observable["observable"],
                    sample["label"],
                    sample.get("comparison_status", ""),
                    str(bool(sample.get("empty", False))).lower(),
                    table_value(sample.get("integral")),
                    table_value(sample.get("integral_error")),
                    table_value(sample.get("mc_max_rel_uncertainty")),
                    ""
                    if sample.get("mc_zero_bin_count") is None
                    else str(sample["mc_zero_bin_count"]),
                    table_value(sample.get("mc_effective_events")),
                    table_value(sample.get("mc_ess_fraction")),
                    table_value(sample.get("data_integral")),
                    table_value(sample.get("data_integral_error")),
                    table_value(sample.get("integral_pull")),
                    table_value(sample.get("chi2")),
                    "" if sample.get("ndf") is None else str(sample["ndf"]),
                    table_value(sample.get("chi2_ndf")),
                    table_value(sample.get("shape_l1")),
                    table_value(sample.get("mc_data_error_ratio")),
                    sample.get("measurement", {}).get("status", "not_tested"),
                    table_value(sample.get("measurement", {}).get("pvalues", {}).get("chi2")),
                    table_value(sample.get("measurement", {}).get("pvalues", {}).get("integral")),
                )
                lines.append("\t".join(values))
    report_path.with_suffix(".txt").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    cprint(f"Saving validation report: {report_path}", "yellow")


# Run the complete iceplot analysis workflow
def main():
    start_time = time.perf_counter()

    args = parse_args()
    if args.validate and args.report is None:
        raise ValueError("--validate requires --report to retain the comparison results")
    print(args)
    plot_output_dir = pathlib.Path(args.output_dir) / args.output
    # Get standalone observables and the visualization scale

    all_obs = None
    if args.analysis is None:
        all_obs = readers.get_observables(args.obs)
        if len(all_obs) == 0:
            raise Exception(__name__ + ".compute: Observables not found!")
        all_obs, unitscale = change_scale(all_obs=all_obs, args=args)
    else:
        unitscale = change_scale(all_obs={}, args=args)[1]
    # Data: Loop over measurement sets

    json_data = None
    data = None
    if args.analysis is not None:
        filename = args.dataset_file

        ### Read out the datasets
        json_data = args.dataset_card

        validate_mc_hepdata_samples(
            assignments=args.mc_hepdata_sample,
            dataset_sets=json_data["sets"],
            input_count=len(args.hepmc3),
        )
        NSETS = len(json_data["sets"])
        data = [None] * NSETS
        name = [None] * NSETS
        set_titles = [None] * NSETS
        set_mc_scales = [None] * NSETS

        # Triplets
        obs = [None] * NSETS
        pid = [None] * NSETS
        cuts = [None] * NSETS
        set_unit_scales = [None] * NSETS

        for j in range(NSETS):
            plot_name = json_data["sets"][j].get("plotname", json_data["sets"][j]["name"])
            name[j] = (plot_name + "/").replace(" ", "_")
            set_titles[j] = json_data["sets"][j].get("title")
            set_mc_scales[j] = float(json_data["sets"][j].get("mc_scale", 1.0))
            ### Get data and observables

            #
            # --> this will return which observables `obs` we actually need to compute
            #
            set_obs, _ = steering.load_observables(
                json_data["sets"][j]["obs"],
                dataset_path=filename,
                cdir=args.cdir,
            )
            unit_overrides = None
            if not args.unit_cli:
                unit_overrides = {
                    histogram["obs"]: histogram["unit"]
                    for histogram in json_data["sets"][j]["hist"]
                    if "unit" in histogram
                }
                if not unit_overrides:
                    unit_overrides = None
            set_obs, set_unit_scales[j] = change_scale(
                all_obs=set_obs,
                args=args,
                units=unit_overrides,
            )
            if json_data["sets"][j].get("data", True):
                hepdata, data_obs = readers.read_hepdata(
                    dataset=json_data["sets"][j],
                    datapath=json_data["datapath"],
                    datatype=json_data["type"],
                    all_obs=set_obs,
                    cdir=args.cdir,
                    reader=json_data.get("reader"),
                    dataset_path=filename,
                )

                if args.datascale is None:
                    scale = 1.0
                elif len(args.datascale) == 1:
                    scale = args.datascale[0]
                else:
                    scale = args.datascale[j]
                scale = apply_unit_scale(scale, set_unit_scales[j])

                data[j] = plots.histhepdata(
                    hepdata=hepdata,
                    obs=data_obs,
                    density=args.density,
                    density_uncertainty=args.density_uncertainty,
                    scale=scale,
                    label=json_data["sets"][j]["name"],
                    hfunc=json_data["plot"]["data_style"],
                    verbose=args.verbose,
                )
                set_histogram_linestyle(
                    [data[j]],
                    json_data["plot"].get("data_linestyle"),
                )
                obs[j] = data_obs if json_data["sets"][j].get("mc", True) else {}
            else:
                obs[j] = select_mc_only_observables(json_data["sets"][j], set_obs)

            ### Collect PID and CUTS for MC

            pid[j] = json_data["sets"][j]["pid"]
            cuts[j] = resolve_cut_reference(
                value=json_data["sets"][j]["cuts"],
                cdir=args.cdir,
                dataset_path=filename,
            )
    # MC (one or several hepmc3 files)

    mclist = [None] * len(args.hepmc3)
    source_diagnostics = [None] * len(args.hepmc3)

    for i in range(len(mclist)):
        hepmc3file = args.hepmc3files[i]
        source_files = args.hepmc3_sources[i]

        print(f"Reading {hepmc3file} ({len(source_files)} file(s))")

        if args.analysis is None:
            cprint("HepDATA card not given, using CLI setup for (obs,pid,cuts)", "magenta")

            NSETS = 1
            name = [""]
            set_titles = [None]
            set_mc_scales = [1.0]

            # Triplets (encapsulate with [], i.e. a single SET)
            obs = [all_obs]
            pid = [args.pid[i]]
            cuts = [args.cut_references[i]]

        # Scaling
        scale = float(args.mcscale[i]) if args.mcscale is not None else 1.0
        if args.analysis is None:
            scale *= unitscale
            scales = effective_mc_scales(scale, set_mc_scales)
        else:
            scales = [
                apply_unit_scale(scale * set_mc_scales[j], set_unit_scales[j])
                for j in range(len(set_mc_scales))
            ]

        assignment = None if args.mc_hepdata_sample is None else args.mc_hepdata_sample[i]
        set_indices = list(
            dataset_indices_for_mc(
                assignment=assignment,
                dataset_sets=json_data["sets"] if json_data is not None else [{"name": value} for value in name],
            )
        )

        if not set_indices:
            mclist[i] = [None] * NSETS
            continue

        # Find out how many events the source contains and whether weights are constant
        cprint(
            f"iceplot: Reading the number of events and event weights in source: {hepmc3file}",
            "yellow",
        )
        weight_summary, file_summaries = scan_hepmc3_source(
            files=source_files,
            maxevents=args.maxevents,
        )
        nevents = weight_summary["nevents"]
        if weight_summary["maxevents_reached"]:
            cprint(f"iceplot: Maximum events limit = {args.maxevents} reached", "red")
        cprint(f"iceplot: Events found in total: {nevents}", "yellow")
        if nevents is None or nevents == 0:
            err_msg = __name__ + f".compute: Error, no events found in source: {hepmc3file}"
            print(err_msg)
            raise Exception(err_msg)

        resolved_xsmode, xsmode_reason = resolve_xsmode(
            requested_xsmode=args.xsmode, weight_summary=weight_summary
        )
        print_xsmode_resolution(
            hepmc3file=hepmc3file,
            requested_xsmode=args.xsmode,
            resolved_xsmode=resolved_xsmode,
            reason=xsmode_reason,
        )
        warn_header_xsmode_weighted_events(
            hepmc3file=hepmc3file, weight_summary=weight_summary, xsmode=args.xsmode
        )
        warn_sample_xsmode_constant_weights(
            hepmc3file=hepmc3file, weight_summary=weight_summary, xsmode=args.xsmode
        )
        warn_sample_xsmode_overflow_events(
            hepmc3file=hepmc3file, weight_summary=weight_summary, xsmode=args.xsmode
        )

        # 2. Chunked loop [parallel processing]
        same_param = {
            "chunk_range": None,
            "hepmc3file": None,
            "obs": [obs[index] for index in set_indices],
            "pid": [pid[index] for index in set_indices],
            "cuts": [cuts[index] for index in set_indices],
            "xsmode": resolved_xsmode,
            "header_xsection": None,
            "scales": [scales[index] for index in set_indices],
            "k": i,
            "label": f"MC [{os.path.basename(hepmc3file)}]"
            if args.mclabel is None
            else args.mclabel[i],
            "density": args.density and not args.stack,
            "density_uncertainty": args.density_uncertainty,
            "verbose": args.verbose,
            "return_diagnostics": True,
        }

        # Create complete, non-overlapping event ranges across all physical files
        paramlist = []
        chunk_sizes = []
        scanned_files = source_files[: len(file_summaries)]
        for filename, file_summary in zip(scanned_files, file_summaries, strict=True):
            file_events = file_summary["nevents"]
            if file_events == 0:
                continue
            event_partitions, _ = build_event_partitions(
                nevents=file_events,
                chunksize=args.chunksize,
                cores=args.cores,
            )
            header_xsection = (
                readers.terminal_header_xsection(file_summary)
                if resolved_xsmode == "header"
                else None
            )
            for chunk_range in event_partitions:
                paramlist.append(
                    {
                        **copy.deepcopy(same_param),
                        "chunk_range": chunk_range,
                        "hepmc3file": filename,
                        "header_xsection": header_xsection,
                    }
                )
                chunk_sizes.append(chunk_range[1] - chunk_range[0] + 1)

        event_workers = min(args.cores, len(paramlist))
        chunk_mode = "automatic" if args.chunksize == "auto" else "fixed"
        cprint(
            f"iceplot: Event chunking [{args.hepmc3[i]}]: {chunk_mode}, "
            f"{len(paramlist)} chunk(s), {min(chunk_sizes)}--{max(chunk_sizes)} events/chunk, "
            f"{event_workers} worker(s)",
            "yellow",
        )

        with multiprocessing.Pool(processes=event_workers) as pool:
            results = pool.map(readers.parallel_wrapper, paramlist)

        # 3. Fuse histograms from each chunk
        diagnostics = aggregate_source_diagnostics(
            chunks=[result["diagnostics"] for result in results],
            xsmode=resolved_xsmode,
        )
        source_diagnostics[i] = {**diagnostics, "after": [None] * NSETS}
        for index, record in zip(set_indices, diagnostics["after"], strict=True):
            source_diagnostics[i]["after"][index] = record
        print_source_summary_table(
            source=hepmc3file,
            summary=diagnostics,
            names=[name[index] for index in set_indices],
            output_dir=plot_output_dir,
            table_name=f"table_source_summary_{i:03d}",
        )
        hists = [result["histograms"] for result in results]
        mclist[i] = [None] * NSETS
        for index, histograms in zip(set_indices, plots.fuse_worker_chunk_outputs(hists), strict=True):
            mclist[i][index] = histograms
        if json_data is not None:
            set_histogram_linestyle(
                mclist[i],
                json_data["plot"].get("mc_linestyle"),
            )
        if not args.stack:
            print_mc_data_comparison_tables(
                source=hepmc3file,
                mc_sets=[mclist[i][index] for index in set_indices],
                data_sets=None if data is None else [data[index] for index in set_indices],
                names=[name[index] for index in set_indices],
                unit=args.unit,
                normalization="unit_density" if args.density else "cross_section",
                output_dir=plot_output_dir,
                table_prefix=f"table_histogram_summary_mc_{i:03d}",
            )

    if len(mclist) == 0 and data is not None:
        print_mc_data_comparison_tables(
            source=None,
            mc_sets=None,
            data_sets=data,
            names=name,
            unit=args.unit,
            normalization="unit_density" if args.density else "cross_section",
            output_dir=plot_output_dir,
            table_prefix="table_histogram_summary_data",
        )
    elif args.stack:
        summed_mc_sets = build_summed_mc_process_sets(
            args=args,
            mclist=mclist,
            names=name,
            dataset_card=json_data,
        )
        print_mc_data_comparison_tables(
            source="MC process stack",
            mc_sets=summed_mc_sets,
            data_sets=data,
            names=name,
            unit=args.unit,
            normalization="unit_density" if args.density else "cross_section",
            output_dir=plot_output_dir,
            table_prefix="table_histogram_summary_stack",
        )

    if args.report is not None:
        report = build_validation_report(
            args=args,
            mclist=mclist,
            data=data,
            names=name,
            dataset_card=json_data,
            source_diagnostics=source_diagnostics,
        )
        write_validation_report(report=report, path=args.report)

    dataset_sets = json_data["sets"] if json_data is not None else [{"name": name[0]}]
    plot_settings = json_data["plot"] if json_data is not None else {}

    # Over each independent set or shared plot group
    for plot_name, set_indices in dataset_plot_groups(dataset_sets):
        ### Plotting

        observable_keys = None
        grouped_mc = []
        grouped_data = []
        grouped_titles = []
        color_index = 0
        data_series_count = (
            0 if data is None else sum(data[index] is not None for index in set_indices)
        )
        data_color = "black" if data_series_count == 1 else None
        for set_index in set_indices:
            source = next((sample[set_index] for sample in mclist if sample[set_index]), None)
            source = source or (data[set_index] if data else None) or {}
            if not source:
                continue
            current_keys = tuple(source.keys())
            if observable_keys is None:
                observable_keys = current_keys
            elif current_keys != observable_keys:
                raise ValueError(f"iceplot: plot group {plot_name!r} has incompatible observables")

            mc_indices = mc_indices_for_dataset(
                assignments=args.mc_hepdata_sample,
                dataset=dataset_sets[set_index] if args.analysis is not None else {},
                input_count=len(mclist),
            )
            mc_sets = [mclist[mc_index][set_index] for mc_index in mc_indices]
            data_set = (
                data[set_index]
                if args.analysis is not None and data[set_index] is not None
                else None
            )
            color_index = append_colored_plot_channel(
                grouped_mc=grouped_mc,
                grouped_data=grouped_data,
                mc_sets=mc_sets,
                data_set=data_set,
                color_index=color_index,
                data_color=data_color,
                match_colors=plot_settings.get("match_colors", False),
            )
            if set_titles[set_index]:
                grouped_titles.append(set_titles[set_index])

        group_title = grouped_titles[0] if len(set(grouped_titles)) == 1 else None

        same_param = {
            "mc": grouped_mc,
            "data": grouped_data,
            "name": plot_name,
            "title": group_title,
            "args": args,
            "ratio_plot": plot_group_ratio_enabled(json_data, set_indices),
            "mc_first": plot_settings.get("mc_first", False),
        }

        # Create parameter dictionaries for each process
        paramlist = [
            dict(copy.deepcopy(same_param), OBS=observable)
            for observable in (observable_keys or ())
        ]

        with multiprocessing.Pool(processes=args.cores) as pool:
            pool.map(create_figures_parallel, paramlist)
    # Save objects out
    # pickle.dump({'mclist': mclist, 'datalist': datalist, 'args': args})
    elapsed_seconds = time.perf_counter() - start_time
    cprint(iceruntime.format_done_message("iceplot", elapsed_seconds), "yellow")
    if args.validate:
        check_physics(report)


if __name__ == "__main__":
    main()
