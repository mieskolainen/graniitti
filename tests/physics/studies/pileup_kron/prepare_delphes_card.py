#!/usr/bin/env python3
# Prepare a Delphes card with a local fixed PU sample
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""Prepare a Delphes card with a local fixed-PU pile-up sample."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path


# Build and validate the command-line parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, help="Input Delphes TCL card")
    parser.add_argument("--output", required=True, help="Output Delphes TCL card")
    parser.add_argument("--pileup-file", required=True, help="Delphes .pileup file")
    parser.add_argument("--mean-pu", type=int, required=True, help="Mean or fixed pile-up")
    parser.add_argument(
        "--fixed-pu",
        action="store_true",
        help="Set PileUpDistribution 2 so MeanPileUp is the exact event count",
    )
    return parser


# Compute true if a line is an active TCL setter for the requested variable
def is_active_setter(line: str, variable: str) -> bool:
    stripped = line.strip()
    return stripped.startswith(f"set {variable} ") and not stripped.startswith("#")


# Replace or add the pile-up settings inside the Delphes PileUpMerger module
def patch_pileup_module(
    lines: list[str], pileup_file: Path, mean_pu: int, fixed_pu: bool
) -> list[str]:
    patched: list[str] = []
    inside_module = False
    saw_distribution = False

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("module PileUpMerger PileUpMerger"):
            inside_module = True
            saw_distribution = False
            patched.append(line)
            continue

        if inside_module and stripped == "}":
            if fixed_pu and not saw_distribution:
                patched.append("  # fixed pile-up multiplicity for deterministic PU overlay\n")
                patched.append("  set PileUpDistribution 2\n")
            inside_module = False
            patched.append(line)
            continue

        if inside_module and is_active_setter(line, "PileUpFile"):
            patched.append(f"  set PileUpFile {pileup_file}\n")
            continue

        if inside_module and is_active_setter(line, "MeanPileUp"):
            patched.append(f"  set MeanPileUp {mean_pu}\n")
            continue

        if inside_module and is_active_setter(line, "PileUpDistribution"):
            saw_distribution = True
            if fixed_pu:
                patched.append("  set PileUpDistribution 2\n")
            else:
                patched.append(line)
            continue

        patched.append(line)

    return patched


# Compute the Delphes branch name from an active TreeWriter line
def tree_writer_branch_name(line: str) -> str | None:
    stripped = line.strip()
    if stripped.startswith("#"):
        return None
    tokens = stripped.split()
    if len(tokens) == 5 and tokens[0] == "add" and tokens[1] == "Branch":
        return tokens[3]
    return None


# Compute TreeWriter branches needed by detector-level graph studies
def graph_input_branch_lines() -> list[str]:
    return [
        "  add Branch HCal/eflowTracks EFlowTrackAll Track\n",
        "  add Branch PhotonEnergySmearing/eflowPhotons EFlowPhoton Tower\n",
        "  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower\n",
    ]


# Insert graph-study EFlow branches into the Delphes TreeWriter module
def patch_tree_writer_graph_branches(lines: list[str]) -> list[str]:
    patched: list[str] = []
    inside_module = False
    active_branches: set[str] = set()
    inserted = False
    insert_before_prefixes = (
        "add Branch Photon",
        "add Branch Electron",
        "add Branch Muon",
        "add Branch JetEnergyScale",
    )

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("module TreeWriter TreeWriter"):
            inside_module = True
            active_branches = set()
            inserted = False
            patched.append(line)
            continue

        if inside_module:
            branch_name = tree_writer_branch_name(line)
            if branch_name is not None:
                active_branches.add(branch_name)

            should_insert = not inserted and any(
                stripped.startswith(prefix) for prefix in insert_before_prefixes
            )
            if should_insert:
                for branch_line in graph_input_branch_lines():
                    name = tree_writer_branch_name(branch_line)
                    if name is not None and name not in active_branches:
                        patched.append(branch_line)
                        active_branches.add(name)
                inserted = True

            if stripped == "}":
                if not inserted:
                    for branch_line in graph_input_branch_lines():
                        name = tree_writer_branch_name(branch_line)
                        if name is not None and name not in active_branches:
                            patched.append(branch_line)
                            active_branches.add(name)
                inside_module = False

        patched.append(line)

    return patched


# Rewrite TCL source includes relative to the generated card directory
def patch_relative_sources(lines: list[str], template_dir: Path, output_dir: Path) -> list[str]:
    patched: list[str] = []
    pattern = re.compile(r"^(\s*)source\s+([^\s#]+)(.*)$")
    for line in lines:
        match = pattern.match(line)
        if not match:
            patched.append(line)
            continue
        indent, source_path, suffix = match.groups()
        path = Path(source_path)
        absolute_path = path.resolve() if path.is_absolute() else (template_dir / path).resolve()
        relative_path = os.path.relpath(absolute_path, output_dir.resolve())
        patched.append(f"{indent}source {relative_path}{suffix}\n")
    return patched


# Read, patch, and write one Delphes card
def prepare_card(
    template: Path, output: Path, pileup_file: Path, mean_pu: int, fixed_pu: bool
) -> None:
    if not template.is_file():
        raise FileNotFoundError(f"Delphes template card not found: {template}")
    if not pileup_file.is_file():
        raise FileNotFoundError(f"Delphes pile-up file not found: {pileup_file}")

    lines = template.read_text(encoding="utf-8").splitlines(keepends=True)
    patched = patch_pileup_module(lines, pileup_file.resolve(), mean_pu, fixed_pu)
    patched = patch_tree_writer_graph_branches(patched)
    output.parent.mkdir(parents=True, exist_ok=True)
    patched = patch_relative_sources(patched, template.resolve().parent, output.parent)
    output.write_text("".join(patched), encoding="utf-8")


# Main entry point for the Delphes card preparation tool
def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    prepare_card(
        template=Path(args.template),
        output=Path(args.output),
        pileup_file=Path(args.pileup_file),
        mean_pu=args.mean_pu,
        fixed_pu=args.fixed_pu,
    )
    print(f"Prepared Delphes card: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
