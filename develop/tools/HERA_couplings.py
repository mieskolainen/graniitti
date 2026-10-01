#!/usr/bin/env python3
# Derive photoproduction parameters from HERA data
#
# The SOFT proton coupling is fixed and HERA data determine the gamma-Pomeron-vector coupling
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import pathlib
import sys

from core.io import files, serialize

# Resolve the development package when this file is executed directly
if __package__ in (None, ""):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from develop.tools import lib
from develop.tools.lib import common
from develop.tools.lib.hera import couplings, fit, model, tune


# Parse command-line options for the HERA derivation
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Fit vector meson photoproduction couplings to HERA data')
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        help='output directory (default: dated figs/HERA directory)',
    )
    parser.add_argument(
        "--jpsi-diagnostics",
        type=pathlib.Path,
        help='write fit inputs, pulls and sensitivity checks to JSON',
    )
    parser.add_argument(
        "--lhcb-fit",
        type=pathlib.Path,
        nargs="?",
        const=lib.SETTINGS,
        help='fit HERA and LHCb using the hera settings',
    )
    parser.add_argument(
        "--write-tune", type=pathlib.Path, help='write a separate tune with fitted parameters'
    )
    parser.add_argument(
        "--g-p-pp",
        type=float,
        default=None,
        help='override g_Ppp [GeV^-1] for reports',
    )
    parser.add_argument(
        "--tune-general",
        type=pathlib.Path,
        default=model.REPOSITORY_ROOT / "modeldata" / "TUNE0" / "GENERAL.json",
        help='GENERAL.json with the SOFT proton vertex',
    )
    parser.add_argument(
        "--channel",
        choices=("all", "heavy", "rho", "phi", "jpsi", "psi2s", "upsilon1s", "upsilon2s", "upsilon3s"),
        default="all",
        help="print one vector-meson channel or all channels",
    )
    parser.add_argument(
        "--format",
        choices=("table", "json", "cards"),
        default="table",
        help='table with update preview, or JSON/cards output',
    )
    parser.add_argument(
        "--push",
        action="store_true",
        help='preview and confirm card updates',
    )
    return parser.parse_args()


# Execute the HERA photoproduction derivation
def main() -> int:
    args = parse_args()
    if (args.push or args.write_tune is not None or args.lhcb_fit is not None) and args.g_p_pp is not None:
        raise ValueError(
            "--push cannot be combined with --g-p-pp because the beam residue is derived from the mapped SOFT exchange"
        )
    proton = model.load_proton(args.tune_general)
    if args.g_p_pp is not None and (not math.isfinite(args.g_p_pp) or args.g_p_pp <= 0.0):
        raise ValueError("--g-p-pp must be finite and positive")
    if args.g_p_pp is not None:
        proton = dataclasses.replace(
            proton, beam_residue_per_gev=args.g_p_pp, source=proton.source + " with --g-p-pp override"
        )

    settings = serialize.load_json_file(args.lhcb_fit or lib.SETTINGS)["hera"]
    config = settings if args.lhcb_fit else None
    selected = []
    states = couplings.channels(proton, config)
    for channel in states:
        if args.channel not in ("all", channel.key) and not (
            args.channel == "heavy" and channel.pdg in (443, 100443, 553, 100553, 200553)
        ):
            continue
        selected.append((channel, model.derive(channel, proton)))

    if args.write_tune is not None:
        tune.write(args.tune_general.parent, args.write_tune, selected)

    from develop.tools.lib.hera import report

    directory = args.output_dir or files.dated_directory(model.REPOSITORY_ROOT / "figs" / "HERA")
    payload = couplings.output(selected, proton)
    report.write(payload, couplings.cards(payload), directory, proton, settings["diagnostics"])
    print(f"HERA diagnostics: {directory.resolve()}", file=sys.stderr)

    if args.jpsi_diagnostics is not None:
        args.jpsi_diagnostics.parent.mkdir(parents=True, exist_ok=True)
        args.jpsi_diagnostics.write_text(
            json.dumps(
                fit.jpsi_diagnostics(proton, next(channel.dlog for channel in states if channel.pdg == 443)), indent=2
            )
            + "\n"
        )

    if args.push or args.format == "table":
        report.print_tables(selected, proton)
        if args.g_p_pp is None and (args.push or (args.lhcb_fit is None and args.write_tune is None)):
            print(f"\nUpdate preview: {args.tune_general.parent}")
            tune.push_cards(args.tune_general, selected)
    elif args.format == "json":
        print(json.dumps(payload, indent=2, ensure_ascii=False))
    else:
        print(common.dumps(couplings.cards(payload)))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, TypeError, ValueError) as exc:
        print(f"HERA_couplings.py: {exc}", file=sys.stderr)
        raise SystemExit(1) from None
