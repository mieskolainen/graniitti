#!/usr/bin/env python3
# Summarize Delphes jet collections against GenJet
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""Summarize Delphes jet collections against GenJet before adding a new reducer."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import uproot
from core.kinematics import fastjet

try:
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False

    # Compute the original function when numba is unavailable
    def njit(*_args: object, **_kwargs: object) -> object:
        if len(_args) == 1 and callable(_args[0]) and not _kwargs:
            return _args[0]

        def decorator(function: object) -> object:
            return function

        return decorator


@dataclass
class Jet:
    pt: float
    eta: float
    phi: float
    mass: float
    index: int


@dataclass
class Match:
    event: int
    collection: str
    gen_pt: float
    gen_eta: float
    gen_phi: float
    reco_pt: float
    reco_eta: float
    reco_phi: float
    dr: float


@dataclass
class Candidate:
    pt: float
    eta: float
    phi: float
    mass: float
    charge: int
    vertex: int
    kind: str
    index: int


@dataclass
class WeightedParticle:
    pt: float
    eta: float
    phi: float
    mass: float


@dataclass
class GraphJetRecord:
    jet: Jet
    hard_fraction: float
    hard_pt: float
    total_pt: float
    log_odds: float = 0.0
    primary_charged_pt: float = 0.0
    pileup_charged_pt: float = 0.0
    neutral_pt: float = 0.0
    neutral_hard_pt: float = 0.0


@dataclass
class GraphJetComponents:
    hard_fraction: float
    hard_pt: float
    total_pt: float
    primary_charged_pt: float
    pileup_charged_pt: float
    neutral_pt: float
    neutral_hard_pt: float


@dataclass
class FisherGraphModel:
    weights: np.ndarray
    bias: float
    sample_log_prior_odds: float
    hard_mean: np.ndarray
    pileup_mean: np.ndarray
    covariance: np.ndarray
    hard_count: int
    pileup_count: int


# Build and validate the command-line parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input Delphes ROOT file")
    parser.add_argument("--output-dir", required=True, help="Directory for summary outputs")
    parser.add_argument("--reference", default="GenJet", help="Reference jet branch")
    parser.add_argument(
        "--collections", nargs="+", default=["Jet", "JetPUPPI"], help="Reco jet branches"
    )
    parser.add_argument(
        "--match-radius", type=float, default=0.30, help="Maximum DeltaR for gen-reco matching"
    )
    parser.add_argument("--gen-pt-min", type=float, default=30.0, help="Minimum reference jet pT")
    parser.add_argument("--reco-pt-min", type=float, default=15.0, help="Minimum reco jet pT")
    parser.add_argument("--eta-max", type=float, default=4.7, help="Maximum absolute jet eta")
    parser.add_argument("--event-start", type=int, default=0, help="First event index to analyze")
    parser.add_argument(
        "--max-events", type=int, default=0, help="Analyze only the first N events when positive"
    )
    parser.add_argument(
        "--event-indices", default="", help="Comma-separated event indices to analyze"
    )
    parser.add_argument(
        "--event-indices-file", default="", help="File containing event indices to analyze"
    )
    parser.add_argument(
        "--include-graph-kron",
        action="store_true",
        help="Build the standalone Graph-Kron collection",
    )
    parser.add_argument(
        "--graph-kron-collection", default="JetGraphKron", help="Standalone collection name"
    )
    parser.add_argument(
        "--graph-kron-base-axis",
        choices=["posterior", "raw"],
        default="raw",
        help="Anti-kT axes used for accepted graph base jets",
    )
    parser.add_argument(
        "--graph-kron-cluster-mode",
        choices=["posterior", "calibrated", "filter", "weighted"],
        default="posterior",
        help="Standalone graph anti-kT output mode",
    )
    parser.add_argument(
        "--graph-kron-candidate-pt-min", type=float, default=0.50, help="Minimum EFlow candidate pT"
    )
    parser.add_argument(
        "--graph-kron-max-candidates",
        type=int,
        default=1200,
        help="Maximum graph candidates per event",
    )
    parser.add_argument(
        "--graph-kron-graph-radius", type=float, default=0.45, help="Candidate graph radius"
    )
    parser.add_argument(
        "--graph-kron-graph-sigma", type=float, default=0.25, help="Candidate graph Gaussian width"
    )
    parser.add_argument(
        "--graph-kron-radiation-radius",
        type=float,
        default=1.20,
        help="Long-range charged-neutral radiation radius",
    )
    parser.add_argument(
        "--graph-kron-radiation-core",
        type=float,
        default=0.04,
        help="Angular core for radiation correlation",
    )
    parser.add_argument(
        "--graph-kron-radiation-angular-power",
        type=float,
        default=1.0,
        help="Angular power in radiation correlation",
    )
    parser.add_argument(
        "--graph-kron-radiation-pt-power",
        type=float,
        default=1.0,
        help="Charged-pT power in radiation correlation",
    )
    parser.add_argument(
        "--graph-kron-radiation-edge-strength",
        type=float,
        default=0.08,
        help="Long-range radiation edge strength",
    )
    parser.add_argument(
        "--graph-kron-radiation-edge-max-neighbours",
        type=int,
        default=16,
        help="Maximum radiation edges per neutral",
    )
    parser.add_argument(
        "--graph-kron-radiation-unary-strength",
        type=float,
        default=2.0,
        help="Radiation unary support strength",
    )
    parser.add_argument(
        "--graph-kron-anchor-strength",
        type=float,
        default=25.0,
        help="Charged-track terminal strength",
    )
    parser.add_argument(
        "--graph-kron-local-support-strength",
        type=float,
        default=1.5,
        help="Neutral local charged-support strength",
    )
    parser.add_argument(
        "--graph-kron-fisher-mode",
        choices=["off", "sample"],
        default="sample",
        help="Fit and apply charged-track Fisher likelihood terms",
    )
    parser.add_argument(
        "--graph-kron-fisher-strength",
        type=float,
        default=2.0,
        help="Neutral Fisher likelihood pseudo-count strength",
    )
    parser.add_argument(
        "--graph-kron-fisher-shrinkage",
        type=float,
        default=0.20,
        help="Fisher covariance diagonal shrinkage",
    )
    parser.add_argument(
        "--graph-kron-fisher-min-class-count",
        type=int,
        default=25,
        help="Minimum charged labels per Fisher class",
    )
    parser.add_argument(
        "--graph-kron-fisher-event-prior-weight",
        type=float,
        default=0.35,
        help="Event prior weight in Fisher odds",
    )
    parser.add_argument(
        "--graph-kron-orphan-neutral-strength",
        type=float,
        default=0.0,
        help="Neutral-orphan hard prior strength",
    )
    parser.add_argument(
        "--graph-kron-orphan-neutral-pt-scale",
        type=float,
        default=15.0,
        help="Neutral pT scale for orphan hard prior",
    )
    parser.add_argument(
        "--graph-kron-orphan-neutral-pileup-scale",
        type=float,
        default=1.0,
        help="Pile-up support scale for orphan hard prior",
    )
    parser.add_argument(
        "--graph-kron-sink-strength",
        type=float,
        default=0.20,
        help="Unsupported-neutral sink strength",
    )
    parser.add_argument(
        "--graph-kron-iterations", type=int, default=20, help="Label propagation iterations"
    )
    parser.add_argument(
        "--graph-kron-pu-terminals", type=int, default=12, help="Maximum pile-up vertex terminals"
    )
    parser.add_argument(
        "--graph-kron-weight-min",
        type=float,
        default=0.03,
        help="Minimum neutral hard-scatter weight",
    )
    parser.add_argument(
        "--graph-kron-jet-radius", type=float, default=0.40, help="Hybrid jet clustering radius"
    )
    parser.add_argument(
        "--graph-kron-cluster-pt-min", type=float, default=5.0, help="Minimum built hybrid jet pT"
    )
    parser.add_argument(
        "--graph-kron-output-scale",
        type=float,
        default=0.45,
        help="Final hybrid jet four-vector scale",
    )
    parser.add_argument(
        "--graph-kron-neutral-pileup-scale",
        type=float,
        default=1.0,
        help="Charged-pileup to neutral-pileup subtraction scale",
    )
    parser.add_argument(
        "--graph-kron-subtracted-hard-strength",
        type=float,
        default=0.0,
        help="Recovery blend strength for charged-subtracted neutral hard pT",
    )
    parser.add_argument(
        "--graph-kron-min-jet-hard-fraction",
        type=float,
        default=0.70,
        help="Minimum graph hard fraction for anti-kT jets",
    )
    parser.add_argument(
        "--graph-kron-min-jet-hard-pt",
        type=float,
        default=18.0,
        help="Minimum graph hard pT for anti-kT jets",
    )
    parser.add_argument(
        "--graph-kron-residual-pu-fraction",
        type=float,
        default=0.35,
        help="Residual non-hard graph component kept in anti-kT jets",
    )
    parser.add_argument(
        "--graph-kron-recovery-min-hard-fraction",
        type=float,
        default=0.20,
        help="Minimum graph hard fraction for recovery jets",
    )
    parser.add_argument(
        "--graph-kron-recovery-min-hard-pt",
        type=float,
        default=15.0,
        help="Minimum graph hard pT for recovery jets",
    )
    parser.add_argument(
        "--graph-kron-recovery-pt-min",
        type=float,
        default=15.0,
        help="Minimum calibrated pT for recovered graph jets",
    )
    parser.add_argument(
        "--graph-kron-recovery-dedupe-radius",
        type=float,
        default=0.28,
        help="Dedupe radius for recovered graph jets",
    )
    parser.add_argument(
        "--graph-kron-recovery-max-add",
        type=int,
        default=1,
        help="Maximum graph recovery jets per event",
    )
    parser.add_argument(
        "--graph-kron-recovery-residual-fraction",
        type=float,
        default=0.35,
        help="Residual non-hard graph component kept in recovery jets",
    )
    parser.add_argument(
        "--graph-kron-recovery-score-pileup-penalty",
        type=float,
        default=0.50,
        help="Pile-up penalty in recovery jet ranking",
    )
    parser.add_argument(
        "--graph-kron-recovery-log-odds-min",
        type=float,
        default=2.0,
        help="Minimum graph likelihood log-odds for recovery jets",
    )
    parser.add_argument(
        "--graph-kron-recovery-likelihood-weight",
        type=float,
        default=1.0,
        help="Recovery jet Fisher likelihood weight",
    )
    parser.add_argument(
        "--graph-kron-recovery-output-scale",
        type=float,
        default=0.45,
        help="Scale applied to recovered graph jets",
    )
    parser.add_argument(
        "--graph-kron-latent-recovery-raw-pt-min",
        type=float,
        default=70.0,
        help="Minimum raw graph jet pT for latent recovery",
    )
    parser.add_argument(
        "--graph-kron-latent-recovery-min-hard-pt",
        type=float,
        default=1.0,
        help="Minimum graph hard support for latent recovery",
    )
    parser.add_argument(
        "--graph-kron-latent-recovery-max-add",
        type=int,
        default=1,
        help="Maximum latent graph recovery jets per event",
    )
    parser.add_argument(
        "--graph-kron-latent-recovery-output-scale",
        type=float,
        default=0.30,
        help="Scale applied to latent raw graph jets",
    )
    parser.add_argument(
        "--no-progress", action="store_true", help="Disable event-level progress bars"
    )
    return parser


# Validate the baseline analysis options before reading the ROOT file
def validate_args(args: argparse.Namespace) -> None:
    if args.match_radius <= 0.0:
        raise ValueError("match-radius must be positive")
    if args.gen_pt_min < 0.0:
        raise ValueError("gen-pt-min must be non-negative")
    if args.reco_pt_min < 0.0:
        raise ValueError("reco-pt-min must be non-negative")
    if args.eta_max <= 0.0:
        raise ValueError("eta-max must be positive")
    if args.event_start < 0:
        raise ValueError("event-start must be non-negative")
    if args.max_events < 0:
        raise ValueError("max-events must be non-negative")
    parse_event_indices(args)
    if args.graph_kron_candidate_pt_min < 0.0:
        raise ValueError("graph-kron-candidate-pt-min must be non-negative")
    if args.graph_kron_max_candidates <= 0:
        raise ValueError("graph-kron-max-candidates must be positive")
    if args.graph_kron_graph_radius <= 0.0:
        raise ValueError("graph-kron-graph-radius must be positive")
    if args.graph_kron_graph_sigma <= 0.0:
        raise ValueError("graph-kron-graph-sigma must be positive")
    if args.graph_kron_radiation_radius <= 0.0:
        raise ValueError("graph-kron-radiation-radius must be positive")
    if args.graph_kron_radiation_core <= 0.0:
        raise ValueError("graph-kron-radiation-core must be positive")
    if args.graph_kron_radiation_angular_power <= 0.0:
        raise ValueError("graph-kron-radiation-angular-power must be positive")
    if args.graph_kron_radiation_pt_power < 0.0:
        raise ValueError("graph-kron-radiation-pt-power must be non-negative")
    if args.graph_kron_radiation_edge_strength < 0.0:
        raise ValueError("graph-kron-radiation-edge-strength must be non-negative")
    if args.graph_kron_radiation_edge_max_neighbours < 0:
        raise ValueError("graph-kron-radiation-edge-max-neighbours must be non-negative")
    if args.graph_kron_radiation_unary_strength < 0.0:
        raise ValueError("graph-kron-radiation-unary-strength must be non-negative")
    if args.graph_kron_anchor_strength <= 0.0:
        raise ValueError("graph-kron-anchor-strength must be positive")
    if args.graph_kron_local_support_strength < 0.0:
        raise ValueError("graph-kron-local-support-strength must be non-negative")
    if args.graph_kron_fisher_strength < 0.0:
        raise ValueError("graph-kron-fisher-strength must be non-negative")
    if not 0.0 <= args.graph_kron_fisher_shrinkage <= 1.0:
        raise ValueError("graph-kron-fisher-shrinkage must be in [0, 1]")
    if args.graph_kron_fisher_min_class_count < 0:
        raise ValueError("graph-kron-fisher-min-class-count must be non-negative")
    if not 0.0 <= args.graph_kron_fisher_event_prior_weight <= 1.0:
        raise ValueError("graph-kron-fisher-event-prior-weight must be in [0, 1]")
    if args.graph_kron_orphan_neutral_strength < 0.0:
        raise ValueError("graph-kron-orphan-neutral-strength must be non-negative")
    if args.graph_kron_orphan_neutral_pt_scale <= 0.0:
        raise ValueError("graph-kron-orphan-neutral-pt-scale must be positive")
    if args.graph_kron_orphan_neutral_pileup_scale < 0.0:
        raise ValueError("graph-kron-orphan-neutral-pileup-scale must be non-negative")
    if args.graph_kron_sink_strength < 0.0:
        raise ValueError("graph-kron-sink-strength must be non-negative")
    if args.graph_kron_iterations <= 0:
        raise ValueError("graph-kron-iterations must be positive")
    if args.graph_kron_pu_terminals < 0:
        raise ValueError("graph-kron-pu-terminals must be non-negative")
    if not 0.0 <= args.graph_kron_weight_min <= 1.0:
        raise ValueError("graph-kron-weight-min must be in [0, 1]")
    if args.graph_kron_jet_radius <= 0.0:
        raise ValueError("graph-kron-jet-radius must be positive")
    if args.graph_kron_cluster_pt_min < 0.0:
        raise ValueError("graph-kron-cluster-pt-min must be non-negative")
    if args.graph_kron_output_scale <= 0.0:
        raise ValueError("graph-kron-output-scale must be positive")
    if args.graph_kron_neutral_pileup_scale < 0.0:
        raise ValueError("graph-kron-neutral-pileup-scale must be non-negative")
    if not 0.0 <= args.graph_kron_subtracted_hard_strength <= 1.0:
        raise ValueError("graph-kron-subtracted-hard-strength must be in [0, 1]")
    if not 0.0 <= args.graph_kron_min_jet_hard_fraction <= 1.0:
        raise ValueError("graph-kron-min-jet-hard-fraction must be in [0, 1]")
    if args.graph_kron_min_jet_hard_pt < 0.0:
        raise ValueError("graph-kron-min-jet-hard-pt must be non-negative")
    if not 0.0 <= args.graph_kron_residual_pu_fraction <= 1.0:
        raise ValueError("graph-kron-residual-pu-fraction must be in [0, 1]")
    if not 0.0 <= args.graph_kron_recovery_min_hard_fraction <= 1.0:
        raise ValueError("graph-kron-recovery-min-hard-fraction must be in [0, 1]")
    if args.graph_kron_recovery_min_hard_pt < 0.0:
        raise ValueError("graph-kron-recovery-min-hard-pt must be non-negative")
    if args.graph_kron_recovery_pt_min < 0.0:
        raise ValueError("graph-kron-recovery-pt-min must be non-negative")
    if args.graph_kron_recovery_dedupe_radius <= 0.0:
        raise ValueError("graph-kron-recovery-dedupe-radius must be positive")
    if args.graph_kron_recovery_max_add < 0:
        raise ValueError("graph-kron-recovery-max-add must be non-negative")
    if not 0.0 <= args.graph_kron_recovery_residual_fraction <= 1.0:
        raise ValueError("graph-kron-recovery-residual-fraction must be in [0, 1]")
    if args.graph_kron_recovery_score_pileup_penalty < 0.0:
        raise ValueError("graph-kron-recovery-score-pileup-penalty must be non-negative")
    if args.graph_kron_recovery_likelihood_weight < 0.0:
        raise ValueError("graph-kron-recovery-likelihood-weight must be non-negative")
    if args.graph_kron_recovery_output_scale <= 0.0:
        raise ValueError("graph-kron-recovery-output-scale must be positive")
    if args.graph_kron_latent_recovery_raw_pt_min < 0.0:
        raise ValueError("graph-kron-latent-recovery-raw-pt-min must be non-negative")
    if args.graph_kron_latent_recovery_min_hard_pt < 0.0:
        raise ValueError("graph-kron-latent-recovery-min-hard-pt must be non-negative")
    if args.graph_kron_latent_recovery_max_add < 0:
        raise ValueError("graph-kron-latent-recovery-max-add must be non-negative")
    if args.graph_kron_latent_recovery_output_scale <= 0.0:
        raise ValueError("graph-kron-latent-recovery-output-scale must be positive")


# Compute the wrapped azimuthal angle difference
def delta_phi(phi_a: float, phi_b: float) -> float:
    return math.atan2(math.sin(phi_a - phi_b), math.cos(phi_a - phi_b))


# Compute an iterable wrapped with tqdm when available and enabled
def progress_iter(iterable: object, total: int, description: str, disabled: bool) -> object:
    if disabled:
        return iterable
    try:
        from tqdm import tqdm
    except ImportError:
        return iterable
    return tqdm(iterable, total=total, desc=description, unit="event")


# Compute DeltaR between two eta-phi coordinates
def delta_r(eta_a: float, phi_a: float, eta_b: float, phi_b: float) -> float:
    return math.hypot(eta_a - eta_b, delta_phi(phi_a, phi_b))


# Decode one ROOT basket as a typed numpy array
def decode_basket_data(
    basket: uproot.models.TBasket.Model_TBasket, dtype: str | np.dtype
) -> np.ndarray:
    return np.frombuffer(basket.data.tobytes(), dtype=np.dtype(dtype))


# Compute the count branch name associated with a split Delphes leaf
def count_branch_name(branch: str) -> str:
    return branch.split(".", 1)[0]


# Load a scalar count branch from raw ROOT baskets
def load_count_branch(tree: uproot.behaviors.TBranch.HasBranches, branch: str) -> list[int]:
    if branch not in tree:
        raise KeyError(f"Missing Delphes count branch: {branch}")
    root_branch = tree[branch]
    counts = [0 for _ in range(root_branch.num_entries)]
    for basket_index in range(root_branch.num_baskets):
        start, stop = root_branch.basket_entry_start_stop(basket_index)
        values = decode_basket_data(root_branch.basket(basket_index), ">i4")
        for entry, value in zip(range(int(start), int(stop)), values, strict=False):
            count = int(value)
            if count < 0:
                raise ValueError(f"{branch} count event {entry} is negative: {count}")
            counts[entry] = count
    return counts


# Compute the big-endian storage dtype for a split Delphes numeric branch
def delphes_branch_dtype(
    root_branch: uproot.behaviors.TBranch.TBranch, expected_kind: str
) -> np.dtype:
    typename = str(getattr(root_branch, "typename", ""))
    dtype_by_typename = {
        "float[]": ("float", np.dtype(">f4")),
        "double[]": ("float", np.dtype(">f8")),
        "int32_t[]": ("int", np.dtype(">i4")),
        "int[]": ("int", np.dtype(">i4")),
    }
    branch_name = str(getattr(root_branch, "name", "<unknown>"))
    if typename not in dtype_by_typename:
        raise TypeError(f"Unsupported Delphes branch type for {branch_name}: {typename}")
    actual_kind, dtype = dtype_by_typename[typename]
    if actual_kind != expected_kind:
        raise TypeError(
            f"Delphes branch {branch_name} has type {typename}, expected {expected_kind} branch"
        )
    return dtype


# Compute basket entry offsets in units of decoded branch values
def basket_item_offsets(
    basket: uproot.models.TBasket.Model_TBasket,
    item_size: int,
    counts: list[int],
    start: int,
    stop: int,
) -> list[int]:
    offsets = getattr(basket, "byte_offsets", None)
    if offsets is not None and len(offsets) == int(stop - start) + 1:
        return [int(offset) // item_size for offset in offsets]
    local_counts = counts[int(start) : int(stop)]
    item_offsets = [0]
    for count in local_counts:
        item_offsets.append(item_offsets[-1] + count)
    return item_offsets


# Validate basket offsets before slicing decoded branch data
def validate_item_offsets(
    branch: str, item_offsets: list[int], value_count: int, start: int, stop: int
) -> None:
    expected_offsets = int(stop - start) + 1
    if len(item_offsets) != expected_offsets:
        raise ValueError(
            f"{branch} basket {start}:{stop} has {len(item_offsets)} offsets, expected {expected_offsets}"
        )
    if item_offsets and item_offsets[0] != 0:
        raise ValueError(f"{branch} basket {start}:{stop} starts at item offset {item_offsets[0]}")
    if any(second < first for first, second in zip(item_offsets, item_offsets[1:], strict=False)):
        raise ValueError(f"{branch} basket {start}:{stop} has decreasing item offsets")
    if item_offsets[-1] > value_count:
        raise ValueError(
            f"{branch} basket {start}:{stop} offsets require {item_offsets[-1]} values, decoded {value_count}"
        )


# Load a split Delphes numeric leaf from raw ROOT baskets
def load_jagged_numeric_branch(
    tree: uproot.behaviors.TBranch.HasBranches,
    branch: str,
    expected_kind: str,
) -> list[list[float]] | list[list[int]]:
    if branch not in tree:
        raise KeyError(f"Missing Delphes branch: {branch}")
    root_branch = tree[branch]
    dtype = delphes_branch_dtype(root_branch, expected_kind)
    counts = load_count_branch(tree, count_branch_name(branch))
    events: list[list[float]] | list[list[int]] = [[] for _ in range(root_branch.num_entries)]
    convert = float if expected_kind == "float" else int
    for basket_index in range(root_branch.num_baskets):
        start, stop = root_branch.basket_entry_start_stop(basket_index)
        basket = root_branch.basket(basket_index)
        values = decode_basket_data(basket, dtype)
        item_offsets = basket_item_offsets(basket, dtype.itemsize, counts, int(start), int(stop))
        validate_item_offsets(branch, item_offsets, len(values), int(start), int(stop))
        for local_entry, entry in enumerate(range(int(start), int(stop))):
            first = item_offsets[local_entry]
            last = item_offsets[local_entry + 1]
            events[entry] = [convert(value) for value in values[first:last]]
    return events


# Load a split Delphes float or double leaf from raw ROOT baskets
def load_jagged_float_branch(
    tree: uproot.behaviors.TBranch.HasBranches, branch: str
) -> list[list[float]]:
    events = load_jagged_numeric_branch(tree, branch, "float")
    return [[float(value) for value in event] for event in events]


# Load a split Delphes integer leaf from raw ROOT baskets
def load_jagged_int_branch(
    tree: uproot.behaviors.TBranch.HasBranches, branch: str
) -> list[list[int]]:
    events = load_jagged_numeric_branch(tree, branch, "int")
    return [[int(value) for value in event] for event in events]


# Compute the first available branch name from a list of aliases
def first_existing_branch(tree: uproot.behaviors.TBranch.HasBranches, names: list[str]) -> str:
    for name in names:
        if name in tree:
            return name
    raise KeyError(f"Missing Delphes branch, tried: {', '.join(names)}")


# Validate that related jagged branches have consistent event and item counts
def validate_equal_event_lengths(context: str, branches: dict[str, list[list[object]]]) -> None:
    if not branches:
        return
    first_name = next(iter(branches))
    event_count = len(branches[first_name])
    for name, events in branches.items():
        if len(events) != event_count:
            raise ValueError(
                f"{context} branch {name} has {len(events)} events, expected {event_count}"
            )
    for event_index in range(event_count):
        lengths = {name: len(events[event_index]) for name, events in branches.items()}
        if len(set(lengths.values())) != 1:
            detail = ", ".join(f"{name}={length}" for name, length in lengths.items())
            raise ValueError(f"{context} event {event_index} branch length mismatch: {detail}")


# Validate that collections refer to the same number of events
def validate_equal_event_counts(context: str, collections: dict[str, list[list[object]]]) -> None:
    if not collections:
        return
    first_name = next(iter(collections))
    event_count = len(collections[first_name])
    for name, events in collections.items():
        if len(events) != event_count:
            raise ValueError(
                f"{context} collection {name} has {len(events)} events, expected {event_count}"
            )


# Validate that all branch values are finite numbers
def validate_finite_branch(branch: str, events: list[list[float]]) -> None:
    for event_index, values in enumerate(events):
        for value_index, value in enumerate(values):
            if not math.isfinite(float(value)):
                raise ValueError(
                    f"{branch} event {event_index} item {value_index} is non-finite: {value}"
                )


# Validate that all branch values are non-negative numbers
def validate_nonnegative_branch(branch: str, events: list[list[float]]) -> None:
    for event_index, values in enumerate(events):
        for value_index, value in enumerate(values):
            if float(value) < 0.0:
                raise ValueError(
                    f"{branch} event {event_index} item {value_index} is negative: {value}"
                )


# Convert Delphes branch arrays into per-event Jet objects
def load_jets(tree: uproot.behaviors.TBranch.HasBranches, branch: str) -> list[list[Jet]]:
    names = [f"{branch}.PT", f"{branch}.Eta", f"{branch}.Phi", f"{branch}.Mass"]
    missing = [name for name in names if name not in tree]
    if missing:
        raise KeyError(f"Missing Delphes jet branch arrays: {', '.join(missing)}")
    pts_by_event = load_jagged_float_branch(tree, names[0])
    etas_by_event = load_jagged_float_branch(tree, names[1])
    phis_by_event = load_jagged_float_branch(tree, names[2])
    masses_by_event = load_jagged_float_branch(tree, names[3])
    validate_equal_event_lengths(
        branch,
        {
            names[0]: pts_by_event,
            names[1]: etas_by_event,
            names[2]: phis_by_event,
            names[3]: masses_by_event,
        },
    )
    for name, values in zip(
        names, [pts_by_event, etas_by_event, phis_by_event, masses_by_event], strict=False
    ):
        validate_finite_branch(name, values)
    validate_nonnegative_branch(names[0], pts_by_event)
    events: list[list[Jet]] = []
    for pts, etas, phis, masses in zip(
        pts_by_event, etas_by_event, phis_by_event, masses_by_event, strict=False
    ):
        jets = [
            Jet(float(pt), float(eta), float(phi), float(mass), index)
            for index, (pt, eta, phi, mass) in enumerate(zip(pts, etas, phis, masses, strict=False))
        ]
        events.append(jets)
    return events


# Build hard-scatter truth anti-kT jets from visible stable particles
def load_hard_truth_jets(
    tree: uproot.behaviors.TBranch.HasBranches,
    radius: float,
    pt_min: float,
    eta_max: float,
) -> list[list[Jet]]:
    required = [
        "Particle.PT",
        "Particle.Eta",
        "Particle.Phi",
        "Particle.Mass",
        "Particle.PID",
        "Particle.Status",
        "Particle.IsPU",
    ]
    missing = [name for name in required if name not in tree]
    if missing:
        raise KeyError(f"Missing hard-truth particle branch arrays: {', '.join(missing)}")
    pts_by_event = load_jagged_float_branch(tree, "Particle.PT")
    etas_by_event = load_jagged_float_branch(tree, "Particle.Eta")
    phis_by_event = load_jagged_float_branch(tree, "Particle.Phi")
    masses_by_event = load_jagged_float_branch(tree, "Particle.Mass")
    pids_by_event = load_jagged_int_branch(tree, "Particle.PID")
    statuses_by_event = load_jagged_int_branch(tree, "Particle.Status")
    is_pu_by_event = load_jagged_int_branch(tree, "Particle.IsPU")
    validate_equal_event_lengths(
        "Particle",
        {
            "Particle.PT": pts_by_event,
            "Particle.Eta": etas_by_event,
            "Particle.Phi": phis_by_event,
            "Particle.Mass": masses_by_event,
            "Particle.PID": pids_by_event,
            "Particle.Status": statuses_by_event,
            "Particle.IsPU": is_pu_by_event,
        },
    )
    for name, values in [
        ("Particle.PT", pts_by_event),
        ("Particle.Eta", etas_by_event),
        ("Particle.Phi", phis_by_event),
        ("Particle.Mass", masses_by_event),
    ]:
        validate_finite_branch(name, values)
    validate_nonnegative_branch("Particle.PT", pts_by_event)
    neutrinos = {12, 14, 16}
    jets_by_event: list[list[Jet]] = []
    for event_values in zip(
        pts_by_event,
        etas_by_event,
        phis_by_event,
        masses_by_event,
        pids_by_event,
        statuses_by_event,
        is_pu_by_event,
        strict=False,
    ):
        particles: list[WeightedParticle] = []
        for particle_values in zip(*event_values, strict=False):
            pt, eta, phi, mass, pid, status, is_pu = particle_values
            if int(status) != 1 or int(is_pu) != 0:
                continue
            if abs(int(pid)) in neutrinos:
                continue
            if pt <= 0.0:
                continue
            particles.append(
                WeightedParticle(float(pt), float(eta), float(phi), max(float(mass), 0.0))
            )
        jets_by_event.append(select_jets(cluster_jets(particles, radius, pt_min), pt_min, eta_max))
    return jets_by_event


# Compute a stable eta-phi tile key for local neighbour searches
def eta_phi_tile(eta: float, phi: float, tile_size: float, phi_tiles: int) -> tuple[int, int]:
    eta_bin = math.floor(eta / tile_size)
    wrapped_phi = (phi + math.pi) % (2.0 * math.pi)
    phi_bin = math.floor(wrapped_phi / tile_size) % phi_tiles
    return eta_bin, phi_bin


# Yield neighbouring eta-phi tile keys with periodic phi wrapping
def neighbouring_tiles(tile: tuple[int, int], phi_tiles: int) -> list[tuple[int, int]]:
    eta_bin, phi_bin = tile
    return [
        (eta_bin + delta_eta, (phi_bin + delta_phi_bin) % phi_tiles)
        for delta_eta in (-1, 0, 1)
        for delta_phi_bin in (-1, 0, 1)
    ]


# Convert a pT eta phi mass tuple to Cartesian four-vector components
def four_vector_from_pt_eta_phi_mass(
    pt: float, eta: float, phi: float, mass: float
) -> tuple[float, float, float, float]:
    px = pt * math.cos(phi)
    py = pt * math.sin(phi)
    pz = pt * math.sinh(eta)
    momentum2 = px * px + py * py + pz * pz
    energy = math.sqrt(max(momentum2 + mass * mass, 0.0))
    return px, py, pz, energy


# Cluster with the shared rapidity based anti-kt implementation and E scheme recombination
def cluster_jets(
    particles: list[WeightedParticle], radius: float, pt_min: float
) -> list[Jet]:
    if not particles:
        return []
    vectors = np.asarray([four_vector_from_pt_eta_phi_mass(p.pt, p.eta, p.phi, p.mass)
                          for p in particles], dtype=np.float64)
    jets, _ = fastjet.cluster_antikt(*vectors.T, radius, pt_min, math.inf)
    return [Jet(fastjet.jet_pt(*p[:2]), fastjet.jet_eta(*p[:3]),
                fastjet.jet_phi(*p[:2]), fastjet.jet_mass(*p), index)
            for index, p in enumerate(jets)]


# Compute the primary vertex index from reconstructed vertex information
def primary_vertex_index(vertex_indices: list[int], vertex_sumpt2: list[float]) -> int:
    if not vertex_indices:
        return 0
    if not vertex_sumpt2:
        return int(vertex_indices[0])
    best = max(range(len(vertex_indices)), key=lambda index: vertex_sumpt2[index])
    return int(vertex_indices[best])


# Compute the selected vertex terminals, primary first
def selected_vertex_terminals(
    vertex_indices: list[int], vertex_sumpt2: list[float], primary_vertex: int, max_pu: int
) -> list[int]:
    scored_vertices = [
        (float(vertex_sumpt2[index]) if index < len(vertex_sumpt2) else 0.0, int(vertex))
        for index, vertex in enumerate(vertex_indices)
        if int(vertex) != primary_vertex
    ]
    scored_vertices.sort(reverse=True)
    return [primary_vertex] + [vertex for _score, vertex in scored_vertices[:max_pu]]


# Load reconstructed vertex arrays needed for graph terminals
def load_vertices(
    tree: uproot.behaviors.TBranch.HasBranches,
) -> tuple[list[list[int]], list[list[float]], list[list[float]]]:
    index_branch = first_existing_branch(tree, ["Vertex.Index"])
    z_branch = first_existing_branch(tree, ["Vertex.Z"])
    sumpt2_branch = first_existing_branch(tree, ["Vertex.SumPT2"])
    vertex_indices = load_jagged_int_branch(tree, index_branch)
    vertex_z = load_jagged_float_branch(tree, z_branch)
    vertex_sumpt2 = load_jagged_float_branch(tree, sumpt2_branch)
    validate_equal_event_lengths(
        "Vertex", {index_branch: vertex_indices, z_branch: vertex_z, sumpt2_branch: vertex_sumpt2}
    )
    validate_finite_branch(z_branch, vertex_z)
    validate_finite_branch(sumpt2_branch, vertex_sumpt2)
    validate_nonnegative_branch(sumpt2_branch, vertex_sumpt2)
    return vertex_indices, vertex_z, vertex_sumpt2


# Assign one track to the nearest reconstructed vertex in z
def assign_track_vertex(
    track_z: float, vertex_indices: list[int], vertex_z: list[float], fallback: int
) -> int:
    if not vertex_indices or not vertex_z or not math.isfinite(track_z):
        return fallback
    best = min(range(len(vertex_indices)), key=lambda index: abs(track_z - vertex_z[index]))
    return int(vertex_indices[best])


# Build EFlow candidates for all events
def load_candidates(
    tree: uproot.behaviors.TBranch.HasBranches,
    pt_min: float,
    eta_max: float,
    vertex_indices_by_event: list[list[int]],
    vertex_z_by_event: list[list[float]],
) -> list[list[Candidate]]:
    track_pt = load_jagged_float_branch(tree, first_existing_branch(tree, ["EFlowTrackAll.PT"]))
    track_eta = load_jagged_float_branch(tree, first_existing_branch(tree, ["EFlowTrackAll.Eta"]))
    track_phi = load_jagged_float_branch(tree, first_existing_branch(tree, ["EFlowTrackAll.Phi"]))
    track_mass = load_jagged_float_branch(tree, first_existing_branch(tree, ["EFlowTrackAll.Mass"]))
    track_z = load_jagged_float_branch(tree, first_existing_branch(tree, ["EFlowTrackAll.Z"]))
    track_charge = load_jagged_int_branch(
        tree, first_existing_branch(tree, ["EFlowTrackAll.Charge"])
    )
    track_vertex = load_jagged_int_branch(
        tree, first_existing_branch(tree, ["EFlowTrackAll.VertexIndex"])
    )
    photon_pt = load_jagged_float_branch(tree, first_existing_branch(tree, ["EFlowPhoton.ET"]))
    photon_eta = load_jagged_float_branch(tree, first_existing_branch(tree, ["EFlowPhoton.Eta"]))
    photon_phi = load_jagged_float_branch(tree, first_existing_branch(tree, ["EFlowPhoton.Phi"]))
    neutral_pt = load_jagged_float_branch(
        tree, first_existing_branch(tree, ["EFlowNeutralHadron.ET"])
    )
    neutral_eta = load_jagged_float_branch(
        tree, first_existing_branch(tree, ["EFlowNeutralHadron.Eta"])
    )
    neutral_phi = load_jagged_float_branch(
        tree, first_existing_branch(tree, ["EFlowNeutralHadron.Phi"])
    )
    validate_equal_event_lengths(
        "EFlowTrackAll",
        {
            "PT": track_pt,
            "Eta": track_eta,
            "Phi": track_phi,
            "Mass": track_mass,
            "Z": track_z,
            "Charge": track_charge,
            "VertexIndex": track_vertex,
        },
    )
    validate_equal_event_lengths(
        "EFlowPhoton", {"ET": photon_pt, "Eta": photon_eta, "Phi": photon_phi}
    )
    validate_equal_event_lengths(
        "EFlowNeutralHadron", {"ET": neutral_pt, "Eta": neutral_eta, "Phi": neutral_phi}
    )
    for name, events in [
        ("EFlowTrackAll.PT", track_pt),
        ("EFlowPhoton.ET", photon_pt),
        ("EFlowNeutralHadron.ET", neutral_pt),
    ]:
        validate_finite_branch(name, events)
        validate_nonnegative_branch(name, events)
    validate_finite_branch("EFlowTrackAll.Z", track_z)
    event_count = len(track_pt)
    candidates_by_event: list[list[Candidate]] = []
    for event_index in range(event_count):
        candidates: list[Candidate] = []
        next_index = 0
        vertex_indices = vertex_indices_by_event[event_index]
        vertex_z = vertex_z_by_event[event_index]
        for values in zip(
            track_pt[event_index],
            track_eta[event_index],
            track_phi[event_index],
            track_mass[event_index],
            track_z[event_index],
            track_charge[event_index],
            track_vertex[event_index],
            strict=False,
        ):
            pt, eta, phi, mass, z, charge, vertex = values
            if pt >= pt_min and abs(eta) <= eta_max and math.isfinite(eta) and math.isfinite(phi):
                assigned_vertex = int(vertex)
                if assigned_vertex < 0:
                    assigned_vertex = assign_track_vertex(
                        float(z), vertex_indices, vertex_z, assigned_vertex
                    )
                candidates.append(
                    Candidate(
                        pt,
                        eta,
                        phi,
                        max(mass, 0.0),
                        int(charge),
                        assigned_vertex,
                        "track",
                        next_index,
                    )
                )
                next_index += 1
        for pt, eta, phi in zip(
            photon_pt[event_index], photon_eta[event_index], photon_phi[event_index], strict=False
        ):
            if pt >= pt_min and abs(eta) <= eta_max and math.isfinite(eta) and math.isfinite(phi):
                candidates.append(Candidate(pt, eta, phi, 0.0, 0, -1, "photon", next_index))
                next_index += 1
        for pt, eta, phi in zip(
            neutral_pt[event_index],
            neutral_eta[event_index],
            neutral_phi[event_index],
            strict=False,
        ):
            if pt >= pt_min and abs(eta) <= eta_max and math.isfinite(eta) and math.isfinite(phi):
                candidates.append(Candidate(pt, eta, phi, 0.0, 0, -1, "neutral", next_index))
                next_index += 1
        candidates_by_event.append(candidates)
    return candidates_by_event


# Compute candidates with compact indices after selection
def reindex_candidates(candidates: list[Candidate]) -> list[Candidate]:
    return [
        Candidate(c.pt, c.eta, c.phi, c.mass, c.charge, c.vertex, c.kind, index)
        for index, c in enumerate(candidates)
    ]


# Keep a bounded candidate set with charged anchors and neutral energy
def limit_candidates(
    candidates: list[Candidate], max_candidates: int, primary_vertex: int
) -> list[Candidate]:
    if len(candidates) <= max_candidates:
        return reindex_candidates(candidates)
    primary_charged = [
        candidate
        for candidate in candidates
        if candidate.charge != 0 and candidate.vertex == primary_vertex
    ]
    neutral = [candidate for candidate in candidates if candidate.charge == 0]
    pileup_charged = [
        candidate
        for candidate in candidates
        if candidate.charge != 0 and candidate.vertex != primary_vertex
    ]
    primary_quota = min(len(primary_charged), max_candidates * 45 // 100)
    neutral_quota = min(len(neutral), max_candidates * 35 // 100)
    pileup_quota = min(len(pileup_charged), max_candidates * 20 // 100)
    selected = (
        sorted(primary_charged, key=lambda candidate: candidate.pt, reverse=True)[:primary_quota]
        + sorted(neutral, key=lambda candidate: candidate.pt, reverse=True)[:neutral_quota]
        + sorted(pileup_charged, key=lambda candidate: candidate.pt, reverse=True)[:pileup_quota]
    )
    selected_ids = {id(candidate) for candidate in selected}
    remaining_budget = max_candidates - len(selected)
    if remaining_budget > 0:
        leftovers = [candidate for candidate in candidates if id(candidate) not in selected_ids]
        selected.extend(
            sorted(leftovers, key=lambda candidate: candidate.pt, reverse=True)[:remaining_budget]
        )
    return reindex_candidates(selected)


# Compute compact numeric arrays for candidate hot loops
def candidate_arrays(
    candidates: list[Candidate],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pt = np.asarray([candidate.pt for candidate in candidates], dtype=np.float64)
    eta = np.asarray([candidate.eta for candidate in candidates], dtype=np.float64)
    phi = np.asarray([candidate.phi for candidate in candidates], dtype=np.float64)
    charge = np.asarray([candidate.charge for candidate in candidates], dtype=np.int64)
    vertex = np.asarray([candidate.vertex for candidate in candidates], dtype=np.int64)
    return pt, eta, phi, charge, vertex


# Compute terminal indices for each candidate vertex
def candidate_terminal_targets(vertices: np.ndarray, terminals: list[int]) -> np.ndarray:
    terminal_targets = np.full(len(vertices), len(terminals), dtype=np.int64)
    for index, vertex in enumerate(vertices):
        for terminal_index, terminal in enumerate(terminals):
            if int(vertex) == int(terminal):
                terminal_targets[index] = terminal_index
                break
    return terminal_targets


# Compute wrapped azimuthal separation in numba hot loops
@njit
def delta_phi_numba(phi_a: float, phi_b: float) -> float:
    return math.atan2(math.sin(phi_a - phi_b), math.cos(phi_a - phi_b))


# Compute squared DeltaR in numba hot loops
@njit
def delta_r2_numba(eta_a: float, phi_a: float, eta_b: float, phi_b: float) -> float:
    dphi = delta_phi_numba(phi_a, phi_b)
    deta = eta_a - eta_b
    return deta * deta + dphi * dphi


# Compute the soft-collinear radiation kernel in numba hot loops
@njit
def radiation_kernel_numba(
    charged_pt: float, dr: float, core: float, angular_power: float, pt_power: float
) -> float:
    pt_weight = max(charged_pt, 0.0) ** pt_power
    angle = math.sqrt(dr * dr + core * core)
    return math.log1p(pt_weight / max(angle**angular_power, 1.0e-12))


# Compute local and radiation terminal supports from labelled charged particles
@njit
def terminal_support_numba(
    pt: np.ndarray,
    eta: np.ndarray,
    phi: np.ndarray,
    charge: np.ndarray,
    terminal_target: np.ndarray,
    terminal_count: int,
    local_radius: float,
    local_sigma: float,
    radiation_radius: float,
    radiation_core: float,
    radiation_angular_power: float,
    radiation_pt_power: float,
) -> tuple[np.ndarray, np.ndarray]:
    size = len(pt)
    local_support = np.zeros((size, terminal_count + 1), dtype=np.float64)
    radiation_support = np.zeros((size, terminal_count + 1), dtype=np.float64)
    local_radius2 = local_radius * local_radius
    radiation_radius2 = radiation_radius * radiation_radius
    sigma2 = local_sigma * local_sigma
    for index in range(size):
        for charged_index in range(size):
            if charge[charged_index] == 0 or charged_index == index:
                continue
            dr2 = delta_r2_numba(eta[index], phi[index], eta[charged_index], phi[charged_index])
            target = terminal_target[charged_index]
            if dr2 <= local_radius2:
                local_support[index, target] += pt[charged_index] * math.exp(-0.5 * dr2 / sigma2)
            if dr2 <= radiation_radius2:
                radiation_support[index, target] += radiation_kernel_numba(
                    pt[charged_index],
                    math.sqrt(dr2),
                    radiation_core,
                    radiation_angular_power,
                    radiation_pt_power,
                )
    return local_support, radiation_support


# Compute a CSR graph with local and bounded radiation edges
@njit
def candidate_graph_csr_numba(
    pt: np.ndarray,
    eta: np.ndarray,
    phi: np.ndarray,
    charge: np.ndarray,
    radius: float,
    sigma: float,
    radiation_radius: float,
    radiation_core: float,
    radiation_angular_power: float,
    radiation_pt_power: float,
    radiation_edge_strength: float,
    radiation_edge_max_neighbours: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    size = len(pt)
    counts = np.zeros(size, dtype=np.int64)
    radius2 = radius * radius
    radiation_radius2 = radiation_radius * radiation_radius
    for index in range(size):
        for other_index in range(index + 1, size):
            if (
                delta_r2_numba(eta[index], phi[index], eta[other_index], phi[other_index])
                <= radius2
            ):
                counts[index] += 1
                counts[other_index] += 1
    if radiation_edge_strength > 0.0:
        for index in range(size):
            if charge[index] != 0:
                continue
            within_count = 0
            for charged_index in range(size):
                if charge[charged_index] == 0:
                    continue
                if (
                    delta_r2_numba(eta[index], phi[index], eta[charged_index], phi[charged_index])
                    <= radiation_radius2
                ):
                    within_count += 1
            if radiation_edge_max_neighbours > 0 and within_count > radiation_edge_max_neighbours:
                within_count = radiation_edge_max_neighbours
            counts[index] += within_count
            for charged_index in range(size):
                if charge[charged_index] == 0:
                    continue
                if (
                    delta_r2_numba(eta[index], phi[index], eta[charged_index], phi[charged_index])
                    <= radiation_radius2
                ):
                    counts[charged_index] += 1 if radiation_edge_max_neighbours == 0 else 0
        if radiation_edge_max_neighbours > 0:
            for index in range(size):
                if charge[index] != 0:
                    continue
                top_weights = np.full(radiation_edge_max_neighbours, -1.0, dtype=np.float64)
                top_indices = np.full(radiation_edge_max_neighbours, -1, dtype=np.int64)
                top_count = 0
                for charged_index in range(size):
                    if charge[charged_index] == 0:
                        continue
                    dr2 = delta_r2_numba(
                        eta[index], phi[index], eta[charged_index], phi[charged_index]
                    )
                    if dr2 > radiation_radius2:
                        continue
                    weight = radiation_edge_strength * radiation_kernel_numba(
                        pt[charged_index],
                        math.sqrt(dr2),
                        radiation_core,
                        radiation_angular_power,
                        radiation_pt_power,
                    )
                    if top_count < radiation_edge_max_neighbours:
                        top_weights[top_count] = weight
                        top_indices[top_count] = charged_index
                        top_count += 1
                    else:
                        min_position = 0
                        min_weight = top_weights[0]
                        for top_index in range(1, radiation_edge_max_neighbours):
                            if top_weights[top_index] < min_weight:
                                min_position = top_index
                                min_weight = top_weights[top_index]
                        if weight > min_weight:
                            top_weights[min_position] = weight
                            top_indices[min_position] = charged_index
                for top_index in range(top_count):
                    counts[top_indices[top_index]] += 1
    indptr = np.zeros(size + 1, dtype=np.int64)
    for index in range(size):
        indptr[index + 1] = indptr[index] + counts[index]
    indices = np.zeros(indptr[size], dtype=np.int64)
    weights = np.zeros(indptr[size], dtype=np.float64)
    cursor = indptr.copy()
    sigma2 = sigma * sigma
    for index in range(size):
        for other_index in range(index + 1, size):
            dr2 = delta_r2_numba(eta[index], phi[index], eta[other_index], phi[other_index])
            if dr2 > radius2:
                continue
            weight = math.exp(-0.5 * dr2 / sigma2)
            position = cursor[index]
            indices[position] = other_index
            weights[position] = weight
            cursor[index] += 1
            position = cursor[other_index]
            indices[position] = index
            weights[position] = weight
            cursor[other_index] += 1
    if radiation_edge_strength > 0.0:
        for index in range(size):
            if charge[index] != 0:
                continue
            if radiation_edge_max_neighbours <= 0:
                for charged_index in range(size):
                    if charge[charged_index] == 0:
                        continue
                    dr2 = delta_r2_numba(
                        eta[index], phi[index], eta[charged_index], phi[charged_index]
                    )
                    if dr2 > radiation_radius2:
                        continue
                    weight = radiation_edge_strength * radiation_kernel_numba(
                        pt[charged_index],
                        math.sqrt(dr2),
                        radiation_core,
                        radiation_angular_power,
                        radiation_pt_power,
                    )
                    position = cursor[index]
                    indices[position] = charged_index
                    weights[position] = weight
                    cursor[index] += 1
                    position = cursor[charged_index]
                    indices[position] = index
                    weights[position] = weight
                    cursor[charged_index] += 1
                continue
            top_weights = np.full(radiation_edge_max_neighbours, -1.0, dtype=np.float64)
            top_indices = np.full(radiation_edge_max_neighbours, -1, dtype=np.int64)
            top_count = 0
            for charged_index in range(size):
                if charge[charged_index] == 0:
                    continue
                dr2 = delta_r2_numba(eta[index], phi[index], eta[charged_index], phi[charged_index])
                if dr2 > radiation_radius2:
                    continue
                weight = radiation_edge_strength * radiation_kernel_numba(
                    pt[charged_index],
                    math.sqrt(dr2),
                    radiation_core,
                    radiation_angular_power,
                    radiation_pt_power,
                )
                if top_count < radiation_edge_max_neighbours:
                    top_weights[top_count] = weight
                    top_indices[top_count] = charged_index
                    top_count += 1
                else:
                    min_position = 0
                    min_weight = top_weights[0]
                    for top_index in range(1, radiation_edge_max_neighbours):
                        if top_weights[top_index] < min_weight:
                            min_position = top_index
                            min_weight = top_weights[top_index]
                    if weight > min_weight:
                        top_weights[min_position] = weight
                        top_indices[min_position] = charged_index
            for top_index in range(top_count):
                charged_index = top_indices[top_index]
                weight = top_weights[top_index]
                position = cursor[index]
                indices[position] = charged_index
                weights[position] = weight
                cursor[index] += 1
                position = cursor[charged_index]
                indices[position] = index
                weights[position] = weight
                cursor[charged_index] += 1
    return indptr, indices, weights


# Solve graph label propagation over a CSR graph
@njit
def solve_label_propagation_csr_numba(
    indptr: np.ndarray,
    indices: np.ndarray,
    edge_weights: np.ndarray,
    unary_support: np.ndarray,
    strengths: np.ndarray,
    iterations: int,
) -> np.ndarray:
    labels = unary_support.copy()
    size = unary_support.shape[0]
    terminals = unary_support.shape[1]
    for _iteration in range(iterations):
        updated = np.empty_like(labels)
        denominators = strengths.copy()
        for index in range(size):
            for terminal in range(terminals):
                updated[index, terminal] = strengths[index] * unary_support[index, terminal]
            for edge_index in range(indptr[index], indptr[index + 1]):
                neighbour = indices[edge_index]
                weight = edge_weights[edge_index]
                denominators[index] += weight
                for terminal in range(terminals):
                    updated[index, terminal] += weight * labels[neighbour, terminal]
        for index in range(size):
            row_sum = 0.0
            denominator = max(denominators[index], 1.0e-12)
            for terminal in range(terminals):
                value = updated[index, terminal] / denominator
                if value < 0.0:
                    value = 0.0
                elif value > 1.0:
                    value = 1.0
                labels[index, terminal] = value
                row_sum += value
            row_sum = max(row_sum, 1.0e-12)
            for terminal in range(terminals):
                labels[index, terminal] /= row_sum
    return labels


# Compute the soft-collinear radiation correlation kernel
def radiation_kernel(
    charged_pt: float,
    dr: float,
    core: float,
    angular_power: float,
    pt_power: float,
) -> float:
    pt_weight = max(charged_pt, 0.0) ** pt_power
    angle = math.sqrt(dr * dr + core * core)
    return math.log1p(pt_weight / max(angle**angular_power, 1.0e-12))


# Add long-range charged-neutral radiation edges to a graph
def add_radiation_edges(
    adjacency: list[list[tuple[int, float]]],
    candidates: list[Candidate],
    radius: float,
    core: float,
    angular_power: float,
    pt_power: float,
    edge_strength: float,
    max_neighbours: int,
) -> None:
    if edge_strength <= 0.0:
        return
    phi_tiles = max(1, int(math.ceil(2.0 * math.pi / radius)))
    charged_tiles: dict[tuple[int, int], list[int]] = {}
    for index, candidate in enumerate(candidates):
        if candidate.charge == 0:
            continue
        key = eta_phi_tile(candidate.eta, candidate.phi, radius, phi_tiles)
        charged_tiles.setdefault(key, []).append(index)
    for neutral_index, neutral in enumerate(candidates):
        if neutral.charge != 0:
            continue
        edges: list[tuple[float, int]] = []
        key = eta_phi_tile(neutral.eta, neutral.phi, radius, phi_tiles)
        for neighbour_key in neighbouring_tiles(key, phi_tiles):
            for charged_index in charged_tiles.get(neighbour_key, []):
                charged = candidates[charged_index]
                dr = delta_r(neutral.eta, neutral.phi, charged.eta, charged.phi)
                if dr > radius:
                    continue
                weight = edge_strength * radiation_kernel(
                    charged.pt, dr, core, angular_power, pt_power
                )
                edges.append((weight, charged_index))
        if max_neighbours > 0 and len(edges) > max_neighbours:
            edges = sorted(edges, reverse=True)[:max_neighbours]
        for weight, charged_index in edges:
            adjacency[neutral_index].append((charged_index, weight))
            adjacency[charged_index].append((neutral_index, weight))


# Build tiled local graph adjacency for candidates
def build_candidate_graph(
    candidates: list[Candidate],
    radius: float,
    sigma: float,
    radiation_radius: float,
    radiation_core: float,
    radiation_angular_power: float,
    radiation_pt_power: float,
    radiation_edge_strength: float,
    radiation_edge_max_neighbours: int,
) -> list[list[tuple[int, float]]]:
    phi_tiles = max(1, int(math.ceil(2.0 * math.pi / radius)))
    tiles: dict[tuple[int, int], list[int]] = {}
    adjacency: list[list[tuple[int, float]]] = [[] for _ in candidates]
    for index, candidate in enumerate(candidates):
        key = eta_phi_tile(candidate.eta, candidate.phi, radius, phi_tiles)
        for neighbour_key in neighbouring_tiles(key, phi_tiles):
            for other_index in tiles.get(neighbour_key, []):
                other = candidates[other_index]
                dr = delta_r(candidate.eta, candidate.phi, other.eta, other.phi)
                if dr > radius:
                    continue
                weight = math.exp(-0.5 * dr * dr / (sigma * sigma))
                adjacency[index].append((other_index, weight))
                adjacency[other_index].append((index, weight))
        tiles.setdefault(key, []).append(index)
    add_radiation_edges(
        adjacency,
        candidates,
        radiation_radius,
        radiation_core,
        radiation_angular_power,
        radiation_pt_power,
        radiation_edge_strength,
        radiation_edge_max_neighbours,
    )
    return adjacency


# Compute charged-terminal support around each candidate
def local_terminal_support(
    candidates: list[Candidate],
    terminals: list[int],
    radius: float,
    sigma: float,
) -> np.ndarray:
    terminal_lookup = {vertex: index for index, vertex in enumerate(terminals)}
    sink_index = len(terminals)
    support = np.zeros((len(candidates), len(terminals) + 1), dtype=float)
    phi_tiles = max(1, int(math.ceil(2.0 * math.pi / radius)))
    charged_tiles: dict[tuple[int, int], list[tuple[int, Candidate]]] = {}
    for index, candidate in enumerate(candidates):
        if candidate.charge == 0:
            continue
        key = eta_phi_tile(candidate.eta, candidate.phi, radius, phi_tiles)
        charged_tiles.setdefault(key, []).append((index, candidate))
    for index, candidate in enumerate(candidates):
        key = eta_phi_tile(candidate.eta, candidate.phi, radius, phi_tiles)
        for neighbour_key in neighbouring_tiles(key, phi_tiles):
            for charged_index, charged_candidate in charged_tiles.get(neighbour_key, []):
                if charged_index == index:
                    continue
                dr = delta_r(
                    candidate.eta, candidate.phi, charged_candidate.eta, charged_candidate.phi
                )
                if dr > radius:
                    continue
                target = terminal_lookup.get(charged_candidate.vertex, sink_index)
                support[index, target] += charged_candidate.pt * math.exp(
                    -0.5 * dr * dr / (sigma * sigma)
                )
    return support


# Compute long-range radiation support around each candidate
def radiation_terminal_support(
    candidates: list[Candidate],
    terminals: list[int],
    radius: float,
    core: float,
    angular_power: float,
    pt_power: float,
) -> np.ndarray:
    terminal_lookup = {vertex: index for index, vertex in enumerate(terminals)}
    sink_index = len(terminals)
    support = np.zeros((len(candidates), len(terminals) + 1), dtype=float)
    phi_tiles = max(1, int(math.ceil(2.0 * math.pi / radius)))
    charged_tiles: dict[tuple[int, int], list[tuple[int, Candidate]]] = {}
    for index, candidate in enumerate(candidates):
        if candidate.charge == 0:
            continue
        key = eta_phi_tile(candidate.eta, candidate.phi, radius, phi_tiles)
        charged_tiles.setdefault(key, []).append((index, candidate))
    for index, candidate in enumerate(candidates):
        key = eta_phi_tile(candidate.eta, candidate.phi, radius, phi_tiles)
        for neighbour_key in neighbouring_tiles(key, phi_tiles):
            for charged_index, charged_candidate in charged_tiles.get(neighbour_key, []):
                if charged_index == index:
                    continue
                dr = delta_r(
                    candidate.eta, candidate.phi, charged_candidate.eta, charged_candidate.phi
                )
                if dr > radius:
                    continue
                target = terminal_lookup.get(charged_candidate.vertex, sink_index)
                support[index, target] += radiation_kernel(
                    charged_candidate.pt, dr, core, angular_power, pt_power
                )
    return support


# Compute terminal supports using the compiled kernel when available
def terminal_support_arrays(
    candidates: list[Candidate],
    terminals: list[int],
    local_radius: float,
    local_sigma: float,
    radiation_radius: float,
    radiation_core: float,
    radiation_angular_power: float,
    radiation_pt_power: float,
) -> tuple[np.ndarray, np.ndarray]:
    if not NUMBA_AVAILABLE:
        return (
            local_terminal_support(candidates, terminals, local_radius, local_sigma),
            radiation_terminal_support(
                candidates,
                terminals,
                radiation_radius,
                radiation_core,
                radiation_angular_power,
                radiation_pt_power,
            ),
        )
    pt, eta, phi, charge, vertex = candidate_arrays(candidates)
    terminal_targets = candidate_terminal_targets(vertex, terminals)
    return terminal_support_numba(
        pt,
        eta,
        phi,
        charge,
        terminal_targets,
        len(terminals),
        local_radius,
        local_sigma,
        radiation_radius,
        radiation_core,
        radiation_angular_power,
        radiation_pt_power,
    )


# Compute a compiled CSR candidate graph for the main event path
def build_candidate_graph_csr(
    candidates: list[Candidate],
    radius: float,
    sigma: float,
    radiation_radius: float,
    radiation_core: float,
    radiation_angular_power: float,
    radiation_pt_power: float,
    radiation_edge_strength: float,
    radiation_edge_max_neighbours: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pt, eta, phi, charge, _vertex = candidate_arrays(candidates)
    return candidate_graph_csr_numba(
        pt,
        eta,
        phi,
        charge,
        radius,
        sigma,
        radiation_radius,
        radiation_core,
        radiation_angular_power,
        radiation_pt_power,
        radiation_edge_strength,
        radiation_edge_max_neighbours,
    )


# Compute the Fisher feature vector for one EFlow candidate
def fisher_feature_vector(
    candidate: Candidate,
    primary_support: float,
    pileup_support: float,
    radiation_primary_support: float = 0.0,
    radiation_pileup_support: float = 0.0,
) -> np.ndarray:
    epsilon = 1.0e-3
    return np.asarray(
        [
            math.log1p(max(candidate.pt, 0.0)),
            math.log1p(max(primary_support, 0.0)),
            math.log1p(max(pileup_support, 0.0)),
            math.log((max(primary_support, 0.0) + epsilon) / (max(pileup_support, 0.0) + epsilon)),
            math.log1p(max(radiation_primary_support, 0.0)),
            math.log1p(max(radiation_pileup_support, 0.0)),
            math.log(
                (max(radiation_primary_support, 0.0) + epsilon)
                / (max(radiation_pileup_support, 0.0) + epsilon)
            ),
            abs(candidate.eta),
        ],
        dtype=float,
    )


# Compute a numerically stable logistic transform
def logistic(value: float) -> float:
    if value >= 50.0:
        return 1.0
    if value <= -50.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-value))


# Compute a finite logit value for a probability
def finite_logit(value: float) -> float:
    clipped = min(max(value, 1.0e-6), 1.0 - 1.0e-6)
    return math.log(clipped / (1.0 - clipped))


# Estimate a regularized Fisher linear discriminant from charged labels
def fit_fisher_model(
    hard_features: list[np.ndarray],
    pileup_features: list[np.ndarray],
    shrinkage: float,
    min_class_count: int,
) -> FisherGraphModel | None:
    if len(hard_features) < min_class_count or len(pileup_features) < min_class_count:
        return None
    hard = np.vstack(hard_features)
    pileup = np.vstack(pileup_features)
    hard_mean = np.mean(hard, axis=0)
    pileup_mean = np.mean(pileup, axis=0)
    hard_centered = hard - hard_mean
    pileup_centered = pileup - pileup_mean
    denominator = max(len(hard_features) + len(pileup_features) - 2, 1)
    covariance = (
        hard_centered.T @ hard_centered + pileup_centered.T @ pileup_centered
    ) / denominator
    diagonal = np.diag(np.diag(covariance))
    covariance = (1.0 - shrinkage) * covariance + shrinkage * diagonal
    covariance = covariance + np.eye(covariance.shape[0], dtype=float) * 1.0e-4
    inverse_covariance = np.linalg.pinv(covariance)
    delta = hard_mean - pileup_mean
    weights = inverse_covariance @ delta
    bias = -0.5 * float(
        hard_mean @ inverse_covariance @ hard_mean - pileup_mean @ inverse_covariance @ pileup_mean
    )
    sample_log_prior_odds = math.log((len(hard_features) + 1.0) / (len(pileup_features) + 1.0))
    return FisherGraphModel(
        weights=weights,
        bias=bias,
        sample_log_prior_odds=sample_log_prior_odds,
        hard_mean=hard_mean,
        pileup_mean=pileup_mean,
        covariance=covariance,
        hard_count=len(hard_features),
        pileup_count=len(pileup_features),
    )


# Compute the event-level charged-track prior odds blended with sample odds
def event_log_prior_odds(
    candidates: list[Candidate],
    primary_vertex: int,
    sample_log_prior_odds: float,
    event_weight: float,
) -> float:
    hard_pt = sum(
        candidate.pt
        for candidate in candidates
        if candidate.charge != 0 and candidate.vertex == primary_vertex
    )
    pileup_pt = sum(
        candidate.pt
        for candidate in candidates
        if candidate.charge != 0 and candidate.vertex != primary_vertex
    )
    event_odds = math.log((hard_pt + 1.0) / (pileup_pt + 1.0))
    return (1.0 - event_weight) * sample_log_prior_odds + event_weight * event_odds


# Compute Fisher hard-scatter log-odds for one feature vector
def fisher_log_odds(model: FisherGraphModel, features: np.ndarray, prior_log_odds: float) -> float:
    return float(model.weights @ features + model.bias + prior_log_odds)


# Estimate sample Fisher parameters from charged EFlow candidates
def estimate_fisher_graph_model(
    candidates_by_event: list[list[Candidate]],
    vertex_indices_by_event: list[list[int]],
    vertex_sumpt2_by_event: list[list[float]],
    args: argparse.Namespace,
) -> FisherGraphModel | None:
    if args.graph_kron_fisher_mode == "off" or args.graph_kron_fisher_strength <= 0.0:
        return None
    hard_features: list[np.ndarray] = []
    pileup_features: list[np.ndarray] = []
    event_iter = zip(
        candidates_by_event, vertex_indices_by_event, vertex_sumpt2_by_event, strict=False
    )
    for candidates, vertex_indices, vertex_sumpt2 in progress_iter(
        event_iter,
        len(candidates_by_event),
        "Graph-Fisher fit",
        args.no_progress,
    ):
        primary = primary_vertex_index(vertex_indices, vertex_sumpt2)
        limited = limit_candidates(candidates, args.graph_kron_max_candidates, primary)
        terminals = selected_vertex_terminals(
            vertex_indices, vertex_sumpt2, primary, args.graph_kron_pu_terminals
        )
        support, radiation_support = terminal_support_arrays(
            limited,
            terminals,
            args.graph_kron_graph_radius,
            args.graph_kron_graph_sigma,
            args.graph_kron_radiation_radius,
            args.graph_kron_radiation_core,
            args.graph_kron_radiation_angular_power,
            args.graph_kron_radiation_pt_power,
        )
        for index, candidate in enumerate(limited):
            if candidate.charge == 0 or candidate.vertex < 0:
                continue
            features = fisher_feature_vector(
                candidate,
                support[index, 0],
                float(np.sum(support[index, 1:])),
                radiation_support[index, 0],
                float(np.sum(radiation_support[index, 1:])),
            )
            if candidate.vertex == primary:
                hard_features.append(features)
            else:
                pileup_features.append(features)
    return fit_fisher_model(
        hard_features,
        pileup_features,
        args.graph_kron_fisher_shrinkage,
        args.graph_kron_fisher_min_class_count,
    )


# Compute candidate-level Fisher log-odds for one event
def candidate_likelihood_log_odds(
    candidates: list[Candidate],
    terminals: list[int],
    fisher_model: FisherGraphModel | None,
    args: argparse.Namespace,
    support_pair: tuple[np.ndarray, np.ndarray] | None = None,
) -> list[float]:
    if fisher_model is None:
        return [0.0 for _candidate in candidates]
    primary = terminals[0]
    if support_pair is None:
        support, radiation_support = terminal_support_arrays(
            candidates,
            terminals,
            args.graph_kron_graph_radius,
            args.graph_kron_graph_sigma,
            args.graph_kron_radiation_radius,
            args.graph_kron_radiation_core,
            args.graph_kron_radiation_angular_power,
            args.graph_kron_radiation_pt_power,
        )
    else:
        support, radiation_support = support_pair
    prior_odds = event_log_prior_odds(
        candidates,
        primary,
        fisher_model.sample_log_prior_odds,
        args.graph_kron_fisher_event_prior_weight,
    )
    log_odds: list[float] = []
    for index, candidate in enumerate(candidates):
        if candidate.charge != 0:
            log_odds.append(8.0 if candidate.vertex == primary else -8.0)
            continue
        features = fisher_feature_vector(
            candidate,
            support[index, 0],
            float(np.sum(support[index, 1:])),
            radiation_support[index, 0],
            float(np.sum(radiation_support[index, 1:])),
        )
        log_odds.append(fisher_log_odds(fisher_model, features, prior_odds))
    return log_odds


# Build local charged-support unary terms for the multi-terminal graph solve
def build_terminal_support(
    candidates: list[Candidate],
    terminals: list[int],
    radius: float,
    sigma: float,
    anchor_strength: float,
    local_support_strength: float,
    sink_strength: float,
    radiation_radius: float = 1.2,
    radiation_core: float = 0.04,
    radiation_angular_power: float = 1.0,
    radiation_pt_power: float = 1.0,
    radiation_unary_strength: float = 0.0,
    fisher_model: FisherGraphModel | None = None,
    fisher_strength: float = 0.0,
    fisher_event_prior_weight: float = 0.0,
    orphan_neutral_strength: float = 0.0,
    orphan_neutral_pt_scale: float = 15.0,
    orphan_neutral_pileup_scale: float = 1.0,
    support_pair: tuple[np.ndarray, np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    terminal_lookup = {vertex: index for index, vertex in enumerate(terminals)}
    sink_index = len(terminals)
    unary_support = np.zeros((len(candidates), len(terminals) + 1), dtype=float)
    strengths = np.full(len(candidates), local_support_strength, dtype=float)
    if support_pair is None:
        local_support, radiation_support = terminal_support_arrays(
            candidates,
            terminals,
            radius,
            sigma,
            radiation_radius,
            radiation_core,
            radiation_angular_power,
            radiation_pt_power,
        )
    else:
        local_support, radiation_support = support_pair
    primary = terminals[0]
    prior_odds = (
        event_log_prior_odds(
            candidates, primary, fisher_model.sample_log_prior_odds, fisher_event_prior_weight
        )
        if fisher_model is not None
        else 0.0
    )
    for index, candidate in enumerate(candidates):
        if candidate.charge != 0:
            target = terminal_lookup.get(candidate.vertex, sink_index)
            unary_support[index, target] = 1.0
            strengths[index] = anchor_strength
            continue
        support = local_support[index].copy() + radiation_unary_strength * radiation_support[index]
        if float(np.sum(support)) <= 0.0:
            support[0] = 0.25
            support[sink_index] = 0.75 + sink_strength
        else:
            support[sink_index] += sink_strength
        local_probability = support / np.sum(support)
        if fisher_model is None or fisher_strength <= 0.0:
            if orphan_neutral_strength > 0.0:
                pileup_evidence = float(
                    np.sum(local_support[index, 1:]) + np.sum(radiation_support[index, 1:])
                )
                orphan_logit = (
                    math.log1p(candidate.pt)
                    - math.log(orphan_neutral_pt_scale)
                    - math.log1p(orphan_neutral_pileup_scale * pileup_evidence)
                )
                orphan_probability = logistic(orphan_neutral_strength * orphan_logit)
                local_probability[0] = max(local_probability[0], orphan_probability)
                local_probability[1:] *= (1.0 - local_probability[0]) / max(
                    float(np.sum(local_probability[1:])), 1.0e-12
                )
            unary_support[index] = local_probability
            continue
        features = fisher_feature_vector(
            candidate,
            local_support[index, 0],
            float(np.sum(local_support[index, 1:])),
            radiation_support[index, 0],
            float(np.sum(radiation_support[index, 1:])),
        )
        hard_probability = logistic(fisher_log_odds(fisher_model, features, prior_odds))
        if orphan_neutral_strength > 0.0:
            pileup_evidence = float(
                np.sum(local_support[index, 1:]) + np.sum(radiation_support[index, 1:])
            )
            orphan_logit = (
                math.log1p(candidate.pt)
                - math.log(orphan_neutral_pt_scale)
                - math.log1p(orphan_neutral_pileup_scale * pileup_evidence)
            )
            hard_probability = max(
                hard_probability, logistic(orphan_neutral_strength * orphan_logit)
            )
        nonhard_probability = local_probability.copy()
        nonhard_probability[0] = 0.0
        nonhard_sum = float(np.sum(nonhard_probability))
        if nonhard_sum <= 0.0:
            nonhard_probability[sink_index] = 1.0
            nonhard_sum = 1.0
        nonhard_probability /= nonhard_sum
        fisher_probability = np.zeros_like(local_probability)
        fisher_probability[0] = hard_probability
        fisher_probability[1:] = (1.0 - hard_probability) * nonhard_probability[1:]
        total_strength = local_support_strength + fisher_strength
        unary_support[index] = (
            local_support_strength * local_probability + fisher_strength * fisher_probability
        ) / total_strength
        strengths[index] = total_strength
    return unary_support, strengths


# Solve the regularized multi-terminal harmonic extension by Jacobi iteration
def solve_label_propagation(
    adjacency: list[list[tuple[int, float]]],
    unary_support: np.ndarray,
    strengths: np.ndarray,
    iterations: int,
) -> np.ndarray:
    labels = unary_support.copy()
    for _iteration in range(iterations):
        updated = strengths[:, None] * unary_support
        denominators = strengths.copy()
        for index, neighbours in enumerate(adjacency):
            for other_index, weight in neighbours:
                updated[index] += weight * labels[other_index]
                denominators[index] += weight
        labels = updated / np.maximum(denominators[:, None], 1.0e-12)
        labels = np.clip(labels, 0.0, 1.0)
        labels /= np.maximum(np.sum(labels, axis=1)[:, None], 1.0e-12)
    return labels


# Build weighted hard-scatter particles from graph labels
def hard_scatter_particles(
    candidates: list[Candidate], labels: np.ndarray, primary_vertex: int, weight_min: float
) -> list[WeightedParticle]:
    particles: list[WeightedParticle] = []
    for candidate, label in zip(candidates, labels, strict=False):
        if candidate.charge != 0:
            if candidate.vertex == primary_vertex:
                particles.append(
                    WeightedParticle(candidate.pt, candidate.eta, candidate.phi, candidate.mass)
                )
            continue
        hard_weight = float(label[0])
        weighted_pt = candidate.pt * hard_weight
        if hard_weight >= weight_min and weighted_pt > 0.0:
            particles.append(
                WeightedParticle(
                    weighted_pt, candidate.eta, candidate.phi, candidate.mass * hard_weight
                )
            )
    return particles


# Build unweighted particles from all graph candidates for anti-kT clustering
def all_candidate_particles(candidates: list[Candidate]) -> list[WeightedParticle]:
    return [
        WeightedParticle(candidate.pt, candidate.eta, candidate.phi, candidate.mass)
        for candidate in candidates
    ]


# Compute per-candidate hard-scatter weights from graph labels
def candidate_hard_weights(
    candidates: list[Candidate], labels: np.ndarray, primary_vertex: int
) -> list[float]:
    weights: list[float] = []
    for candidate, label in zip(candidates, labels, strict=False):
        if candidate.charge != 0:
            weights.append(1.0 if candidate.vertex == primary_vertex else 0.0)
        else:
            weights.append(float(label[0]))
    return weights


# Compute graph hard and charged-neutral component support around one anti-kT jet
def jet_graph_components(
    jet: Jet, candidates: list[Candidate], hard_weights: list[float], radius: float
) -> GraphJetComponents:
    total_pt = 0.0
    hard_pt = 0.0
    primary_charged_pt = 0.0
    pileup_charged_pt = 0.0
    neutral_pt = 0.0
    neutral_hard_pt = 0.0
    for candidate, weight in zip(candidates, hard_weights, strict=False):
        if delta_r(jet.eta, jet.phi, candidate.eta, candidate.phi) > radius:
            continue
        total_pt += candidate.pt
        clipped_weight = max(0.0, min(1.0, weight))
        hard_pt += candidate.pt * clipped_weight
        if candidate.charge != 0:
            if clipped_weight > 0.5:
                primary_charged_pt += candidate.pt
            else:
                pileup_charged_pt += candidate.pt
        else:
            neutral_pt += candidate.pt
            neutral_hard_pt += candidate.pt * clipped_weight
    fraction = hard_pt / total_pt if total_pt > 0.0 else 0.0
    return GraphJetComponents(
        hard_fraction=fraction,
        hard_pt=hard_pt,
        total_pt=total_pt,
        primary_charged_pt=primary_charged_pt,
        pileup_charged_pt=pileup_charged_pt,
        neutral_pt=neutral_pt,
        neutral_hard_pt=neutral_hard_pt,
    )


# Compute graph hard support around one anti-kT jet
def jet_graph_support(
    jet: Jet, candidates: list[Candidate], hard_weights: list[float], radius: float
) -> tuple[float, float, float]:
    components = jet_graph_components(jet, candidates, hard_weights, radius)
    return components.hard_fraction, components.hard_pt, components.total_pt


# Keep and calibrate Python anti-kT jets with graph hard-scatter support
def filter_graph_supported_jets(
    jets: list[Jet],
    candidates: list[Candidate],
    hard_weights: list[float],
    radius: float,
    min_fraction: float,
    min_hard_pt: float,
    residual_fraction: float,
) -> list[Jet]:
    kept_jets: list[Jet] = []
    for jet in jets:
        hard_fraction, hard_pt, _total_pt = jet_graph_support(jet, candidates, hard_weights, radius)
        if hard_fraction >= min_fraction or hard_pt >= min_hard_pt:
            graph_scale = hard_fraction + residual_fraction * (1.0 - hard_fraction)
            kept_jets.append(
                Jet(jet.pt * graph_scale, jet.eta, jet.phi, jet.mass * graph_scale, len(kept_jets))
            )
    return kept_jets


# Calibrate graph-weighted anti-kT jets with local all-candidate support
def calibrate_weighted_graph_jets(
    jets: list[Jet],
    candidates: list[Candidate],
    hard_weights: list[float],
    radius: float,
    min_fraction: float,
    min_hard_pt: float,
    residual_fraction: float,
) -> list[Jet]:
    calibrated_jets: list[Jet] = []
    for jet in jets:
        hard_fraction, hard_pt, total_pt = jet_graph_support(jet, candidates, hard_weights, radius)
        if hard_fraction < min_fraction and hard_pt < min_hard_pt:
            continue
        target_pt = hard_pt + residual_fraction * max(total_pt - hard_pt, 0.0)
        if target_pt <= 0.0 or jet.pt <= 0.0:
            continue
        scale = target_pt / jet.pt
        calibrated_jets.append(
            Jet(target_pt, jet.eta, jet.phi, jet.mass * scale, len(calibrated_jets))
        )
    return calibrated_jets


# Compute an average candidate likelihood log-odds around one jet
def jet_likelihood_log_odds(
    jet: Jet, candidates: list[Candidate], candidate_log_odds: list[float], radius: float
) -> float:
    weighted_sum = 0.0
    total_weight = 0.0
    for candidate, log_odds in zip(candidates, candidate_log_odds, strict=False):
        if delta_r(jet.eta, jet.phi, candidate.eta, candidate.phi) > radius:
            continue
        weight = math.sqrt(max(candidate.pt, 0.0))
        weighted_sum += weight * log_odds
        total_weight += weight
    if total_weight <= 0.0:
        return 0.0
    return weighted_sum / total_weight


# Attach graph hard-scatter support to anti-kT jets
def graph_jet_records(
    jets: list[Jet],
    candidates: list[Candidate],
    hard_weights: list[float],
    radius: float,
    candidate_log_odds: list[float] | None = None,
) -> list[GraphJetRecord]:
    records: list[GraphJetRecord] = []
    log_odds_values = (
        candidate_log_odds if candidate_log_odds is not None else [0.0 for _candidate in candidates]
    )
    for jet in jets:
        components = jet_graph_components(jet, candidates, hard_weights, radius)
        log_odds = jet_likelihood_log_odds(jet, candidates, log_odds_values, radius)
        records.append(
            GraphJetRecord(
                jet=jet,
                hard_fraction=components.hard_fraction,
                hard_pt=components.hard_pt,
                total_pt=components.total_pt,
                log_odds=log_odds,
                primary_charged_pt=components.primary_charged_pt,
                pileup_charged_pt=components.pileup_charged_pt,
                neutral_pt=components.neutral_pt,
                neutral_hard_pt=components.neutral_hard_pt,
            )
        )
    return records


# Compute charged-subtracted neutral hard pT for one graph jet record
def record_subtracted_hard_pt(record: GraphJetRecord, neutral_pileup_scale: float) -> float:
    neutral_hard_pt = max(record.neutral_pt - neutral_pileup_scale * record.pileup_charged_pt, 0.0)
    return max(record.primary_charged_pt + neutral_hard_pt, 0.0)


# Compute the effective hard pT after blending graph and charged-subtracted estimates
def record_effective_hard_pt(
    record: GraphJetRecord, neutral_pileup_scale: float, subtracted_strength: float
) -> float:
    subtracted_hard_pt = record_subtracted_hard_pt(record, neutral_pileup_scale)
    strength = min(max(subtracted_strength, 0.0), 1.0)
    return (1.0 - strength) * record.hard_pt + strength * subtracted_hard_pt


# Compute the effective hard fraction after charged-subtracted blending
def record_effective_hard_fraction(
    record: GraphJetRecord, neutral_pileup_scale: float, subtracted_strength: float
) -> float:
    if record.total_pt <= 0.0:
        return 0.0
    return (
        record_effective_hard_pt(record, neutral_pileup_scale, subtracted_strength)
        / record.total_pt
    )


# Convert one graph support record into a calibrated base jet
def graph_base_jet_from_record(
    record: GraphJetRecord,
    residual_fraction: float,
    output_scale: float,
    neutral_pileup_scale: float,
    subtracted_strength: float,
    index: int,
) -> Jet | None:
    hard_pt = record_effective_hard_pt(record, neutral_pileup_scale, subtracted_strength)
    target_pt = hard_pt + residual_fraction * max(record.total_pt - hard_pt, 0.0)
    if target_pt <= 0.0 or record.jet.pt <= 0.0:
        return None
    mass_scale = output_scale * target_pt / record.jet.pt
    return Jet(
        output_scale * target_pt,
        record.jet.eta,
        record.jet.phi,
        record.jet.mass * mass_scale,
        index,
    )


# Build selected graph posterior base jets from weighted anti-kT records
def graph_base_jets_from_records(
    records: list[GraphJetRecord],
    min_fraction: float,
    min_hard_pt: float,
    residual_fraction: float,
    output_scale: float,
    neutral_pileup_scale: float = 0.0,
    subtracted_strength: float = 0.0,
) -> list[Jet]:
    jets: list[Jet] = []
    for record in records:
        hard_pt = record_effective_hard_pt(record, neutral_pileup_scale, subtracted_strength)
        hard_fraction = hard_pt / record.total_pt if record.total_pt > 0.0 else 0.0
        if hard_fraction < min_fraction and hard_pt < min_hard_pt:
            continue
        jet = graph_base_jet_from_record(
            record,
            residual_fraction,
            output_scale,
            neutral_pileup_scale,
            subtracted_strength,
            len(jets),
        )
        if jet is not None:
            jets.append(jet)
    return jets


# Compute the likelihood-ratio recovery score for one candidate
def recovery_record_log_odds(
    record: GraphJetRecord,
    pileup_penalty: float,
    likelihood_weight: float,
    neutral_pileup_scale: float = 0.0,
    subtracted_strength: float = 0.0,
) -> float:
    hard_pt = record_effective_hard_pt(record, neutral_pileup_scale, subtracted_strength)
    pileup_pt = max(record.total_pt - hard_pt, 0.0)
    hard_fraction = hard_pt / record.total_pt if record.total_pt > 0.0 else 0.0
    fraction_odds = finite_logit(hard_fraction)
    pt_odds = math.log((hard_pt + 1.0) / (pileup_pt + 1.0))
    return (
        fraction_odds
        + pt_odds
        - pileup_penalty * math.log1p(pileup_pt)
        + likelihood_weight * record.log_odds
    )


# Convert one raw graph record into a recovered jet
def recovered_jet_from_record(
    record: GraphJetRecord, args: argparse.Namespace, index: int
) -> Jet | None:
    hard_pt = record_effective_hard_pt(
        record,
        args.graph_kron_neutral_pileup_scale,
        args.graph_kron_subtracted_hard_strength,
    )
    graph_pt = hard_pt + args.graph_kron_recovery_residual_fraction * max(
        record.total_pt - hard_pt, 0.0
    )
    if graph_pt <= 0.0 or record.jet.pt <= 0.0:
        return None
    pt = args.graph_kron_recovery_output_scale * graph_pt
    if pt <= 0.0:
        return None
    mass_scale = pt / record.jet.pt
    return Jet(pt, record.jet.eta, record.jet.phi, record.jet.mass * mass_scale, index)


# Append standalone graph recovery jets from raw anti-kT records
def append_graph_recovery_jets(
    base_jets: list[Jet],
    raw_records: list[GraphJetRecord],
    args: argparse.Namespace,
) -> list[Jet]:
    merged = [Jet(jet.pt, jet.eta, jet.phi, jet.mass, index) for index, jet in enumerate(base_jets)]
    candidates: list[tuple[float, GraphJetRecord]] = []
    for record in raw_records:
        hard_pt = record_effective_hard_pt(
            record,
            args.graph_kron_neutral_pileup_scale,
            args.graph_kron_subtracted_hard_strength,
        )
        hard_fraction = hard_pt / record.total_pt if record.total_pt > 0.0 else 0.0
        if (
            hard_fraction < args.graph_kron_recovery_min_hard_fraction
            and hard_pt < args.graph_kron_recovery_min_hard_pt
        ):
            continue
        if any(
            delta_r(record.jet.eta, record.jet.phi, jet.eta, jet.phi)
            < args.graph_kron_recovery_dedupe_radius
            for jet in merged
        ):
            continue
        log_odds = recovery_record_log_odds(
            record,
            args.graph_kron_recovery_score_pileup_penalty,
            args.graph_kron_recovery_likelihood_weight,
            args.graph_kron_neutral_pileup_scale,
            args.graph_kron_subtracted_hard_strength,
        )
        if log_odds < args.graph_kron_recovery_log_odds_min:
            continue
        candidates.append((log_odds, record))
    candidates.sort(key=lambda item: item[0], reverse=True)
    for _score, record in candidates:
        if len(merged) - len(base_jets) >= args.graph_kron_recovery_max_add:
            break
        recovered = recovered_jet_from_record(record, args, len(merged))
        if recovered is None or recovered.pt < args.graph_kron_recovery_pt_min:
            continue
        if any(
            delta_r(recovered.eta, recovered.phi, jet.eta, jet.phi)
            < args.graph_kron_recovery_dedupe_radius
            for jet in merged
        ):
            continue
        merged.append(Jet(recovered.pt, recovered.eta, recovered.phi, recovered.mass, len(merged)))
    return [Jet(jet.pt, jet.eta, jet.phi, jet.mass, index) for index, jet in enumerate(merged)]


# Append high-raw-pT latent graph jets that have weak but nonzero hard support
def append_latent_graph_recovery_jets(
    base_jets: list[Jet],
    raw_records: list[GraphJetRecord],
    args: argparse.Namespace,
) -> list[Jet]:
    if (
        args.graph_kron_latent_recovery_max_add <= 0
        or args.graph_kron_latent_recovery_raw_pt_min <= 0.0
    ):
        return base_jets
    merged = [Jet(jet.pt, jet.eta, jet.phi, jet.mass, index) for index, jet in enumerate(base_jets)]
    added = 0
    for record in sorted(raw_records, key=lambda item: item.jet.pt, reverse=True):
        if record.jet.pt < args.graph_kron_latent_recovery_raw_pt_min:
            continue
        if record.hard_pt < args.graph_kron_latent_recovery_min_hard_pt:
            continue
        if any(
            delta_r(record.jet.eta, record.jet.phi, jet.eta, jet.phi)
            < args.graph_kron_recovery_dedupe_radius
            for jet in merged
        ):
            continue
        pt = args.graph_kron_latent_recovery_output_scale * record.jet.pt
        if pt < args.reco_pt_min:
            continue
        mass_scale = pt / record.jet.pt if record.jet.pt > 0.0 else 0.0
        merged.append(
            Jet(pt, record.jet.eta, record.jet.phi, record.jet.mass * mass_scale, len(merged))
        )
        added += 1
        if added >= args.graph_kron_latent_recovery_max_add:
            break
    return [Jet(jet.pt, jet.eta, jet.phi, jet.mass, index) for index, jet in enumerate(merged)]


# Apply a common calibration scale to reconstructed jets
def scale_jets(jets: list[Jet], scale: float) -> list[Jet]:
    return [Jet(jet.pt * scale, jet.eta, jet.phi, jet.mass * scale, jet.index) for jet in jets]


# Build one event of Graph-Kron hybrid jets
def build_graph_kron_event(
    candidates: list[Candidate],
    vertex_indices: list[int],
    vertex_sumpt2: list[float],
    args: argparse.Namespace,
    fisher_model: FisherGraphModel | None = None,
) -> list[Jet]:
    primary = primary_vertex_index(vertex_indices, vertex_sumpt2)
    candidates = limit_candidates(candidates, args.graph_kron_max_candidates, primary)
    if not candidates:
        return []
    terminals = selected_vertex_terminals(
        vertex_indices, vertex_sumpt2, primary, args.graph_kron_pu_terminals
    )
    support_pair = terminal_support_arrays(
        candidates,
        terminals,
        args.graph_kron_graph_radius,
        args.graph_kron_graph_sigma,
        args.graph_kron_radiation_radius,
        args.graph_kron_radiation_core,
        args.graph_kron_radiation_angular_power,
        args.graph_kron_radiation_pt_power,
    )
    unary_support, strengths = build_terminal_support(
        candidates,
        terminals,
        args.graph_kron_graph_radius,
        args.graph_kron_graph_sigma,
        args.graph_kron_anchor_strength,
        args.graph_kron_local_support_strength,
        args.graph_kron_sink_strength,
        args.graph_kron_radiation_radius,
        args.graph_kron_radiation_core,
        args.graph_kron_radiation_angular_power,
        args.graph_kron_radiation_pt_power,
        args.graph_kron_radiation_unary_strength,
        fisher_model,
        args.graph_kron_fisher_strength,
        args.graph_kron_fisher_event_prior_weight,
        args.graph_kron_orphan_neutral_strength,
        args.graph_kron_orphan_neutral_pt_scale,
        args.graph_kron_orphan_neutral_pileup_scale,
        support_pair,
    )
    if NUMBA_AVAILABLE:
        indptr, indices, edge_weights = build_candidate_graph_csr(
            candidates,
            args.graph_kron_graph_radius,
            args.graph_kron_graph_sigma,
            args.graph_kron_radiation_radius,
            args.graph_kron_radiation_core,
            args.graph_kron_radiation_angular_power,
            args.graph_kron_radiation_pt_power,
            args.graph_kron_radiation_edge_strength,
            args.graph_kron_radiation_edge_max_neighbours,
        )
        labels = solve_label_propagation_csr_numba(
            indptr, indices, edge_weights, unary_support, strengths, args.graph_kron_iterations
        )
    else:
        adjacency = build_candidate_graph(
            candidates,
            args.graph_kron_graph_radius,
            args.graph_kron_graph_sigma,
            args.graph_kron_radiation_radius,
            args.graph_kron_radiation_core,
            args.graph_kron_radiation_angular_power,
            args.graph_kron_radiation_pt_power,
            args.graph_kron_radiation_edge_strength,
            args.graph_kron_radiation_edge_max_neighbours,
        )
        labels = solve_label_propagation(
            adjacency, unary_support, strengths, args.graph_kron_iterations
        )
    hard_weights = candidate_hard_weights(candidates, labels, primary)
    likelihood_log_odds = candidate_likelihood_log_odds(
        candidates, terminals, fisher_model, args, support_pair
    )
    if args.graph_kron_cluster_mode == "posterior":
        weighted_particles = hard_scatter_particles(
            candidates, labels, primary, args.graph_kron_weight_min
        )
        weighted_jets = cluster_jets(
            weighted_particles, args.graph_kron_jet_radius, args.graph_kron_cluster_pt_min
        )
        weighted_records = graph_jet_records(
            weighted_jets, candidates, hard_weights, args.graph_kron_jet_radius, likelihood_log_odds
        )
        raw_particles = all_candidate_particles(candidates)
        raw_jets = cluster_jets(
            raw_particles, args.graph_kron_jet_radius, args.graph_kron_cluster_pt_min
        )
        raw_records = graph_jet_records(
            raw_jets, candidates, hard_weights, args.graph_kron_jet_radius, likelihood_log_odds
        )
        base_records = raw_records if args.graph_kron_base_axis == "raw" else weighted_records
        base_jets = graph_base_jets_from_records(
            base_records,
            args.graph_kron_min_jet_hard_fraction,
            args.graph_kron_min_jet_hard_pt,
            args.graph_kron_residual_pu_fraction,
            args.graph_kron_output_scale,
            args.graph_kron_neutral_pileup_scale,
            args.graph_kron_subtracted_hard_strength,
        )
        recovered_jets = append_graph_recovery_jets(base_jets, raw_records, args)
        return append_latent_graph_recovery_jets(recovered_jets, raw_records, args)
    if args.graph_kron_cluster_mode == "calibrated":
        particles = hard_scatter_particles(candidates, labels, primary, args.graph_kron_weight_min)
        jets = cluster_jets(
            particles, args.graph_kron_jet_radius, args.graph_kron_cluster_pt_min
        )
        calibrated_jets = calibrate_weighted_graph_jets(
            jets,
            candidates,
            hard_weights,
            args.graph_kron_jet_radius,
            args.graph_kron_min_jet_hard_fraction,
            args.graph_kron_min_jet_hard_pt,
            args.graph_kron_residual_pu_fraction,
        )
        return scale_jets(calibrated_jets, args.graph_kron_output_scale)
    if args.graph_kron_cluster_mode == "filter":
        particles = all_candidate_particles(candidates)
        jets = cluster_jets(
            particles, args.graph_kron_jet_radius, args.graph_kron_cluster_pt_min
        )
        filtered_jets = filter_graph_supported_jets(
            jets,
            candidates,
            hard_weights,
            args.graph_kron_jet_radius,
            args.graph_kron_min_jet_hard_fraction,
            args.graph_kron_min_jet_hard_pt,
            args.graph_kron_residual_pu_fraction,
        )
        return scale_jets(filtered_jets, args.graph_kron_output_scale)
    particles = hard_scatter_particles(candidates, labels, primary, args.graph_kron_weight_min)
    jets = cluster_jets(
        particles, args.graph_kron_jet_radius, args.graph_kron_cluster_pt_min
    )
    return scale_jets(jets, args.graph_kron_output_scale)


# Build Graph-Kron hybrid jets for all events
def build_graph_kron_jets_by_event(
    tree: uproot.behaviors.TBranch.HasBranches,
    args: argparse.Namespace,
) -> list[list[Jet]]:
    vertex_indices_by_event, vertex_z_by_event, vertex_sumpt2_by_event = load_vertices(tree)
    candidates_by_event = load_candidates(
        tree,
        args.graph_kron_candidate_pt_min,
        args.eta_max,
        vertex_indices_by_event,
        vertex_z_by_event,
    )
    event_indices = parse_event_indices(args)
    vertex_indices_by_event = limit_events(
        vertex_indices_by_event, args.max_events, args.event_start, event_indices
    )
    vertex_z_by_event = limit_events(
        vertex_z_by_event, args.max_events, args.event_start, event_indices
    )
    vertex_sumpt2_by_event = limit_events(
        vertex_sumpt2_by_event, args.max_events, args.event_start, event_indices
    )
    candidates_by_event = limit_events(
        candidates_by_event, args.max_events, args.event_start, event_indices
    )
    validate_equal_event_counts(
        "GraphKronInputs",
        {
            "candidates": candidates_by_event,
            "vertices": vertex_indices_by_event,
            "vertex_sumpt2": vertex_sumpt2_by_event,
        },
    )
    fisher_model = estimate_fisher_graph_model(
        candidates_by_event, vertex_indices_by_event, vertex_sumpt2_by_event, args
    )
    jets_by_event: list[list[Jet]] = []
    event_iter = zip(
        candidates_by_event, vertex_indices_by_event, vertex_sumpt2_by_event, strict=False
    )
    for candidates, vertex_indices, vertex_sumpt2 in progress_iter(
        event_iter,
        len(candidates_by_event),
        "Graph-Kron",
        args.no_progress,
    ):
        jets_by_event.append(
            build_graph_kron_event(candidates, vertex_indices, vertex_sumpt2, args, fisher_model)
        )
    return jets_by_event


# Compute jets that pass the kinematic selection
def select_jets(jets: list[Jet], pt_min: float, eta_max: float) -> list[Jet]:
    return [jet for jet in jets if jet.pt >= pt_min and abs(jet.eta) <= eta_max]


# Compute selected event indices from inline text or a file
def parse_event_indices(args: argparse.Namespace) -> list[int] | None:
    inline = str(getattr(args, "event_indices", "") or "").strip()
    file_name = str(getattr(args, "event_indices_file", "") or "").strip()
    if inline and file_name:
        raise ValueError("Use only one of event-indices and event-indices-file")
    text = inline
    if file_name:
        text = Path(file_name).read_text(encoding="utf-8")
    if not text.strip():
        return None
    tokens = text.replace(",", " ").split()
    indices = [int(token) for token in tokens]
    if any(index < 0 for index in indices):
        raise ValueError("event indices must be non-negative")
    if len(set(indices)) != len(indices):
        raise ValueError("event indices must be unique")
    return indices


# Compute the requested deterministic event slice
def limit_events(
    events: list, max_events: int, event_start: int = 0, event_indices: list[int] | None = None
) -> list:
    if event_indices is not None:
        selected = []
        for index in event_indices:
            if index >= len(events):
                raise ValueError(f"event index {index} is outside available range 0:{len(events)}")
            selected.append(events[index])
        return selected
    if event_start > 0:
        events = events[event_start:]
    if max_events <= 0:
        return events
    return events[:max_events]


# Greedily match selected reco jets to selected reference jets
def match_event(
    event_index: int,
    collection: str,
    gen_jets: list[Jet],
    reco_jets: list[Jet],
    match_radius: float,
) -> tuple[list[Match], int]:
    matches: list[Match] = []
    used_reco: set[int] = set()
    for gen in sorted(gen_jets, key=lambda jet: jet.pt, reverse=True):
        best_reco = None
        best_dr = match_radius
        for reco in reco_jets:
            if reco.index in used_reco:
                continue
            dr = delta_r(gen.eta, gen.phi, reco.eta, reco.phi)
            if dr < best_dr:
                best_reco = reco
                best_dr = dr
        if best_reco is None:
            continue
        used_reco.add(best_reco.index)
        matches.append(
            Match(
                event=event_index,
                collection=collection,
                gen_pt=gen.pt,
                gen_eta=gen.eta,
                gen_phi=gen.phi,
                reco_pt=best_reco.pt,
                reco_eta=best_reco.eta,
                reco_phi=best_reco.phi,
                dr=best_dr,
            )
        )
    return matches, len(reco_jets) - len(used_reco)


# Compute a safe mean value from a list
def safe_mean(values: list[float]) -> float:
    if not values:
        return float("nan")
    return float(np.mean(values))


# Count matched jets per event for one collection
def collection_match_counts(matches: list[Match], event_count: int) -> list[int]:
    counts = [0 for _ in range(event_count)]
    for match in matches:
        counts[match.event] += 1
    return counts


# Summarize matches and event counts for one reco collection
def summarize_collection(
    collection: str,
    matches: list[Match],
    gen_counts: list[int],
    reco_counts: list[int],
    fake_counts: list[int],
) -> dict[str, float | int | str]:
    responses = [match.reco_pt / match.gen_pt for match in matches if match.gen_pt > 0.0]
    abs_errors = [abs(response - 1.0) for response in responses]
    dr_values = [match.dr for match in matches]
    fake_fractions = [
        fake / reco if reco > 0 else 0.0
        for fake, reco in zip(fake_counts, reco_counts, strict=False)
    ]
    match_counts = collection_match_counts(matches, len(gen_counts))
    match_efficiencies = [
        min(match_count, gen_count) / gen_count
        for match_count, gen_count in zip(match_counts, gen_counts, strict=False)
        if gen_count > 0
    ]
    total_gen_jets = int(sum(gen_counts))
    return {
        "collection": collection,
        "events": len(gen_counts),
        "gen_jets": total_gen_jets,
        "reco_jets": int(sum(reco_counts)),
        "matched_jets": len(matches),
        "match_efficiency_global": len(matches) / total_gen_jets
        if total_gen_jets > 0
        else float("nan"),
        "match_efficiency_mean": safe_mean(match_efficiencies),
        "fake_fraction_mean": safe_mean(fake_fractions),
        "response_mean": safe_mean(responses),
        "response_median": float(np.median(responses)) if responses else float("nan"),
        "response_width": float(np.std(responses)) if responses else float("nan"),
        "absolute_response_error_mean": safe_mean(abs_errors),
        "match_delta_r_mean": safe_mean(dr_values),
        "reco_jets_per_event_mean": safe_mean([float(count) for count in reco_counts]),
    }


# Write a list of dictionaries as CSV
def write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


# Format one value for markdown table output
def format_value(value: float | int | str) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if not math.isfinite(value):
        return "nan"
    return f"{value:.4g}"


# Write a compact markdown metrics table
def write_markdown(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    headers = list(rows[0].keys())
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(format_value(row[header]) for header in headers) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# Write one row per matched jet
def write_match_rows(path: Path, matches_by_collection: dict[str, list[Match]]) -> None:
    rows: list[dict[str, float | int | str]] = []
    for collection, matches in matches_by_collection.items():
        for match in matches:
            rows.append(
                {
                    "event": match.event,
                    "collection": collection,
                    "gen_pt": match.gen_pt,
                    "gen_eta": match.gen_eta,
                    "gen_phi": match.gen_phi,
                    "reco_pt": match.reco_pt,
                    "reco_eta": match.reco_eta,
                    "reco_phi": match.reco_phi,
                    "delta_r": match.dr,
                    "response": match.reco_pt / match.gen_pt
                    if match.gen_pt > 0.0
                    else float("nan"),
                }
            )
    write_csv(path, rows)


# Analyze one Delphes ROOT file and write baseline metrics
def analyze(args: argparse.Namespace) -> None:
    validate_args(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with uproot.open(
        Path(args.input), handler=uproot.source.file.MultithreadedFileSource
    ) as root_file:
        tree = root_file["Delphes"]
        event_indices = parse_event_indices(args)
        if args.reference == "GenHardJet":
            gen_source = load_hard_truth_jets(
                tree, 0.40, args.graph_kron_cluster_pt_min, args.eta_max
            )
        else:
            gen_source = load_jets(tree, args.reference)
        gen_events = limit_events(gen_source, args.max_events, args.event_start, event_indices)
        reco_events_by_collection = {
            collection: limit_events(
                load_jets(tree, collection), args.max_events, args.event_start, event_indices
            )
            for collection in args.collections
        }
        if args.include_graph_kron:
            reco_events_by_collection[args.graph_kron_collection] = build_graph_kron_jets_by_event(
                tree, args
            )
    validate_equal_event_counts(
        "RecoCollections", {"reference": gen_events, **reco_events_by_collection}
    )
    selected_gen_events = [select_jets(jets, args.gen_pt_min, args.eta_max) for jets in gen_events]
    gen_counts = [len(jets) for jets in selected_gen_events]
    rows: list[dict[str, float | int | str]] = []
    matches_by_collection: dict[str, list[Match]] = {}
    for collection, reco_events in reco_events_by_collection.items():
        selected_reco_events = [
            select_jets(jets, args.reco_pt_min, args.eta_max) for jets in reco_events
        ]
        reco_counts = [len(jets) for jets in selected_reco_events]
        fake_counts: list[int] = []
        matches: list[Match] = []
        for event_index, (gen_jets, reco_jets) in enumerate(
            zip(selected_gen_events, selected_reco_events, strict=False)
        ):
            event_matches, fake_count = match_event(
                event_index, collection, gen_jets, reco_jets, args.match_radius
            )
            matches.extend(event_matches)
            fake_counts.append(fake_count)
        matches_by_collection[collection] = matches
        rows.append(summarize_collection(collection, matches, gen_counts, reco_counts, fake_counts))
    write_csv(output_dir / "delphes_metrics_summary.csv", rows)
    write_markdown(output_dir / "delphes_metrics_summary.md", rows)
    write_match_rows(output_dir / "delphes_matched_jets.csv", matches_by_collection)


# Main entry point for Delphes ROOT baseline analysis
def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    analyze(args)
    print(f"Delphes baseline metrics: {Path(args.output_dir) / 'delphes_metrics_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
