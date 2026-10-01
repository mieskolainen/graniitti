#!/usr/bin/env python3
# Build XP and GP resonance card couplings for two exchange fusion
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
#

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

from core.io.serialize import load_json_file

try:
    import pyjson5 as json5
except ImportError as exc:
    raise SystemExit("resonance_card_builder.py requires pyjson5") from exc

# Resolve the development package when this file is executed directly
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from develop.tools.lib import common, pole, soft_exchange


# Load a JSON-with-comments file
def load_jsonc(path: Path) -> Any:
    return load_json_file(path, loader=json5.load)


# Load particle metadata from built-ins, PDG_EXTRA and RES cards
def load_particles(tune_dir: Path, model: str = "GP") -> dict[int, pole.Particle]:
    particles = pole.particles()
    extra_path = tune_dir / "PDG_EXTRA.json"
    if extra_path.exists():
        for row in load_jsonc(extra_path)["PARAM_PDG"].values():
            pdg = int(row["PDG"])
            particles[pdg] = pole.Particle(
                pdg,
                str(row.get("name", pdg)),
                int(row["spinX2"]),
                int(row["P"]),
                int(row["C"]),
                float(row.get("mass", 0.0)),
                float(row.get("width", 0.0)),
                int(row.get("chargeX3", 0)),
                int(row.get("isospinX2", 0)),
                int(row.get("G", 0)),
            )

    res_dir = tune_dir / "RES"
    if res_dir.exists():
        for path in sorted(res_dir.glob("*.json")):
            row = load_jsonc(path)["PARAM_RES"]
            pdg = int(row["PDG"])
            particles[pdg] = pole.Particle(
                pdg,
                str(row.get("name", path.stem)),
                int(row["spinX2"]),
                int(row["P"]),
                int(row.get("C", 0)),
                float(row["MODELS"][model]["mass"]),
                float(row["MODELS"][model]["width"]),
            )
    return particles


# Load Regge numerics, aliases and signature signs from the tune cards
def load_regge(tune_dir: Path) -> tuple[int, pole.ReggeTable]:
    general = load_jsonc(tune_dir / "GENERAL.json")
    configuration = soft_exchange.load(general)
    numerics = load_jsonc(tune_dir / "NUMERICS.json")
    return int(numerics["NUMERICS_REGGE"]["MMAX"]), pole.ReggeTable(
        groups=[list(row.pdg) for row in configuration.regge_exchanges],
        tau=[
            int(configuration.definitions[row.soft_exchange]["tau"])
            for row in configuration.regge_exchanges
        ],
        pole_spin=[row.pole_spin for row in configuration.regge_exchanges],
    )


# Compute absolute pole LS couplings with one selected active row
def ls_card(
    rows: list[tuple[int, float]], active: tuple[int, float], magnitude: float
) -> list[list[Any]]:
    return [
        [
            l_value,
            pole.clean_number(s_value),
            magnitude if l_value == active[0] and abs(s_value - active[1]) < 1e-9 else 0.0,
            0.0,
        ]
        for l_value, s_value in rows
    ]


# Compute the selected active LS row, defaulting to the first allowed row
def active_row(
    rows: list[tuple[int, float]], active_ls: list[str] | None
) -> tuple[int, float]:
    if not rows:
        raise ValueError("no allowed LS rows")
    if active_ls is None:
        return rows[0]
    active = (int(active_ls[0]), float(active_ls[1]))
    if not any(
        l_value == active[0] and abs(s_value - active[1]) < 1e-9 for l_value, s_value in rows
    ):
        raise ValueError(f"requested active LS row {active} is not allowed")
    return active


# Compute the coupling basis selected for one production model
def read_basis(args: argparse.Namespace) -> str:
    if args.basis not in (None, "g_ls", "helicity"):
        raise ValueError('--basis must be "g_ls" or "helicity"')
    return args.basis or "g_ls"


# Check that one GP exchange is a photon or a null spin analytic alias
def validate_gp_exchange(particle: pole.Particle, regge: pole.ReggeTable) -> None:
    if particle.pdg == 22:
        return
    if particle.spin2 != -1:
        raise ValueError(f"GP exchange PDG {particle.pdg} must have spinX2 = -1")
    if not any(particle.pdg in group for group in regge.groups):
        raise ValueError(f"GP exchange PDG {particle.pdg} is not listed in PARAM_REGGE.EXCHANGES")


# Compute the signature of one mapped analytic exchange
def gp_signature(particle: pole.Particle, regge: pole.ReggeTable) -> int:
    for index, group in enumerate(regge.groups):
        if particle.pdg in group:
            return regge.tau[index]
    raise ValueError(f"GP exchange PDG {particle.pdg} has no trajectory signature")


# Resolve the physical pole with the same charge, isospin and discrete symmetries
def gp_pole(particle: pole.Particle, particles: dict[int, pole.Particle], regge: pole.ReggeTable) -> pole.Particle:
    if particle.pdg == 22:
        return particle
    for index, group in enumerate(regge.groups):
        if particle.pdg not in group:
            continue
        spin = regge.pole_spin[index]
        if regge.tau[index] != (1 if spin % 2 == 0 else -1):
            raise ValueError("GP pole spin is incompatible with the trajectory signature")
        quantum_numbers = (particle.charge3, particle.isospin2, particle.parity,
                           particle.cparity, particle.gparity)
        matches = [
            candidate for pdg in group if (candidate := particles.get(pdg)) is not None
            and candidate.spin2 == 2 * spin
            and (candidate.charge3, candidate.isospin2, candidate.parity,
                 candidate.cparity, candidate.gparity) == quantum_numbers
        ]
        if len(matches) != 1:
            raise ValueError(f"GP exchange PDG {particle.pdg} requires one matching pole alias")
        return matches[0]
    raise ValueError(f"GP exchange PDG {particle.pdg} has no pole mapping")


# Compute the GP Jacob-Wick parity phase for one fusion channel
def gp_parity_phase(
    mother: pole.Particle,
    p1: pole.Particle,
    p2: pole.Particle,
    regge: pole.ReggeTable,
) -> int:
    parity = mother.parity * (1 if mother.spin2 // 2 % 2 == 0 else -1)
    for particle in (p1, p2):
        parity *= particle.parity
        if particle.pdg == 22:
            parity *= 1 if particle.spin2 // 2 % 2 == 0 else -1
        else:
            parity *= gp_signature(particle, regge)
    return parity


# Compute the complete GP pole LS table with the selected active coupling
def gp_ls_rows(
    args: argparse.Namespace, mother: pole.Particle, p1: pole.Particle, p2: pole.Particle
) -> list[list[int | float | None]]:
    if args.active_ls is None:
        raise ValueError("GP g_ls requires --active-ls L S")
    rows = pole.ls_rows(mother, p1, p2, args.c_symmetry, args.p_symmetry)
    active = active_row(rows, args.active_ls)
    magnitude = None if tuple(args.fuse) == (22, 22) else args.channel_mag
    return production_rows(
        ls_card(rows, active, 1.0), magnitude, args.channel_phase, keep_zeros=True
    )


# Compute GP parity and identical exchange orbit with its relative signs
def gp_helicity_orbit(
    m1: int, m2: int, parity: int, exchange: int | None
) -> dict[tuple[int, int], int]:
    partners = [(m1, m2, 1), (-m1, -m2, parity)]
    if exchange is not None:
        partners.extend([(m2, m1, exchange), (-m2, -m1, parity * exchange)])
    orbit: dict[tuple[int, int], int] = {}
    for first, second, phase in partners:
        coordinate = (first, second)
        if coordinate in orbit and orbit[coordinate] != phase:
            raise ValueError("GP parity or identical exchange forces this helicity to zero")
        orbit[coordinate] = phase
    return orbit


# Compute independent GP helicity rows with symmetry partners left implicit
def gp_helicities(
    args: argparse.Namespace,
    mother: pole.Particle,
    p1: pole.Particle,
    p2: pole.Particle,
    mmax: int,
    regge: pole.ReggeTable,
    poles: tuple[pole.Particle, pole.Particle],
) -> list[list[int | float | None]]:
    if not args.helicity_row:
        raise ValueError("GP helicity requires at least one --helicity-row")
    rows: list[list[int | float | None]] = []
    seen_orbits: set[tuple[int, int]] = set()
    parity = gp_parity_phase(mother, p1, p2, regge)
    exchange = None
    if p1.pdg == p2.pdg:
        exponent = (poles[0].spin2 + poles[1].spin2 - mother.spin2) // 2
        exchange = 1 if exponent % 2 == 0 else -1
    for values in args.helicity_row:
        m1_float, m2_float, magnitude, phase = values
        if not pole.is_integer(m1_float) or not pole.is_integer(m2_float):
            raise ValueError("GP helicity m1 and m2 must be integers")
        m1 = int(round(m1_float))
        m2 = int(round(m2_float))
        if abs(m1) > mmax or abs(m2) > mmax:
            raise ValueError("GP helicity row lies outside MMAX")
        if 2 * abs(m1) > poles[0].spin2 or 2 * abs(m2) > poles[1].spin2:
            raise ValueError("GP helicity row exceeds the exchange pole spin")
        if abs(m1 - m2) > 0.5 * mother.spin2:
            raise ValueError("GP helicity row violates |m1-m2| <= JX")
        if (p1.pdg == 22 and abs(m1) != 1) or (p2.pdg == 22 and abs(m2) != 1):
            raise ValueError("GP photon helicity must be transverse")
        orbit = min(gp_helicity_orbit(m1, m2, parity, exchange))
        if orbit in seen_orbits:
            raise ValueError("GP helicity rows repeat one symmetry orbit")
        if magnitude <= 0.0 or not math.isfinite(magnitude) or not math.isfinite(phase):
            raise ValueError("GP helicity magnitude must be finite and positive")
        seen_orbits.add(orbit)
        rows.append([m1, m2, magnitude, common.phase(phase)])
    if tuple(args.fuse) == (22, 22):
        if not math.isclose(rows[0][2], 1.0, rel_tol=0.0, abs_tol=1.0e-12) or abs(rows[0][3]) > 1.0e-12:
            raise ValueError("width-derived gamma gamma raw input needs a unit real first row")
        rows[0][2] = None
    return rows


# Convert a unit shape table into direct production couplings
def production_rows(
    rows: list[list[int | float]], magnitude: float | None, phase: float, *, keep_zeros: bool = False
) -> list[list[int | float | None]]:
    if magnitude is None:
        active = [row for row in rows if float(row[2]) > 0.0]
        if len(active) != 1:
            raise ValueError("width-derived production requires one active row")
        return [
            [row[0], row[1], None if float(row[2]) > 0.0 else 0.0,
             common.phase(float(row[3]) + phase) if float(row[2]) > 0.0 else 0.0]
            for row in (rows if keep_zeros else active)
        ]
    return [
        [
            row[0],
            row[1],
            float(row[2]) * magnitude,
            common.phase(float(row[3]) + phase) if float(row[2]) * magnitude > 0.0 else 0.0,
        ]
        for row in rows
    ]


# Compute XP covariant LS card fragment
def xp_card(
    args: argparse.Namespace,
    mother: pole.Particle,
    p1: pole.Particle,
    p2: pole.Particle,
    rows: list[tuple[int, float]],
    active: tuple[int, float],
    basis: str,
) -> dict[str, Any]:
    diphoton = p1.pdg == 22 and p2.pdg == 22
    magnitude = None if diphoton else args.channel_mag
    block: dict[str, Any] = {
        "basis": basis,
        "Lambda": args.lambda_scale,
        "CP": [args.c_symmetry, args.p_symmetry],
    }
    if basis == "g_ls":
        block["g_ls"] = production_rows(
            ls_card(rows, active, 1.0), magnitude, args.channel_phase, keep_zeros=True
        )
    else:
        block["helicity"] = production_rows(
            pole.helicity_rows(
                mother, p1, p2, rows, active, args.tol, p_symmetry=args.p_symmetry
            ),
            magnitude,
            args.channel_phase,
        )
    key = f"[{p1.pdg},{p2.pdg}]"
    return {"XP": {key: block}}


# Compute the card fragment for the requested model
def build(
    args: argparse.Namespace, particles: dict[int, pole.Particle], mmax: int, regge: pole.ReggeTable
) -> dict[str, Any]:
    mother = pole.Particle(
        args.pdg, args.name, args.spinX2, args.parity, args.cparity, args.mass, args.width
    )
    try:
        p1 = particles[args.fuse[0]]
        p2 = particles[args.fuse[1]]
    except KeyError as exc:
        raise ValueError(f"missing particle metadata for fused PDG {exc.args[0]}") from exc

    basis = read_basis(args)
    if args.model == "XP":
        rows = pole.ls_rows(mother, p1, p2, args.c_symmetry, args.p_symmetry)
        active = active_row(rows, args.active_ls)
        return xp_card(args, mother, p1, p2, rows, active, basis)

    validate_gp_exchange(p1, regge)
    validate_gp_exchange(p2, regge)
    poles = (gp_pole(p1, particles, regge), gp_pole(p2, particles, regge))
    if (
        args.c_symmetry
        and pole.has_cparity(mother)
        and pole.has_cparity(p1)
        and pole.has_cparity(p2)
        and p1.cparity * p2.cparity != mother.cparity
    ):
        raise ValueError("GP trajectory pair violates resonance C parity")

    key = f"[{p1.pdg},{p2.pdg}]"
    block: dict[str, Any] = {
        "basis": basis,
        "CP": [args.c_symmetry, args.p_symmetry],
    }
    if basis == "g_ls":
        block["Lambda"] = args.lambda_scale
        block["g_ls"] = gp_ls_rows(args, mother, *poles)
    else:
        if args.active_ls is not None:
            raise ValueError("GP helicity does not accept --active-ls")
        if args.channel_mag != 1.0 or args.channel_phase != 0.0:
            raise ValueError("GP helicity rows carry their own magnitude and phase")
        block["helicity"] = gp_helicities(args, mother, p1, p2, mmax, regge, poles)
    return {args.model: {key: block}}


# Print generated fragment as summary, g_ls and helicity tables
def print_table(
    args: argparse.Namespace,
    fragment: dict[str, Any],
    particles: dict[int, pole.Particle],
    mmax: int,
) -> None:
    model = fragment[args.model]
    key = next(iter(model))
    channel = [int(value) for value in key[1:-1].split(",")]
    vertex = model[key]
    p1 = particles[channel[0]]
    p2 = particles[channel[1]]
    spin = 0.5 * args.spinX2
    parity_label = "+" if args.parity > 0 else "-"
    cparity_label = "+" if args.cparity > 0 else "-" if args.cparity < 0 else "?"
    common.print_table(
        "Resonance production-card builder",
        ["model", "resonance", "J^PC", "exchange 1", "exchange 2", "MMAX", "basis"],
        [
            [
                args.model,
                args.name,
                f"{pole.clean_number(spin)}^{{{parity_label}{cparity_label}}}",
                f"{p1.name} ({p1.pdg})",
                f"{p2.name} ({p2.pdg})",
                mmax if args.model == "GP" else "n/a",
                vertex["basis"],
            ]
        ],
        right_align={5},
    )

    if vertex["basis"] == "helicity":
        helicity_rows = [
            [
                row[0],
                row[1],
                "width derived" if row[2] is None else f"{float(row[2]):.12g}",
                f"{float(row[3]):.12g}",
            ]
            for row in vertex["helicity"]
        ]
        common.print_table(
            "\nHelicity couplings", ["lambda1", "lambda2", "magnitude", "phase [rad]"], helicity_rows, right_align={0, 1, 2, 3}
        )
        return
    if vertex["basis"] == "g_ls":
        g_ls_rows = [
            [
                row[0],
                row[1],
                "width derived" if row[2] is None else f"{float(row[2]):.12g}",
                f"{float(row[3]):.12g}",
            ]
            for row in vertex["g_ls"]
        ]
        common.print_table("\nCovariant g_ls rows", ["L", "S", "magnitude", "phase [rad]"], g_ls_rows, right_align={0, 1, 2, 3})
        return


# Validate finite CLI values before constructing angular-momentum rows
def validate_args(args: argparse.Namespace, mmax: int) -> None:
    if args.spinX2 < 0:
        raise ValueError("--spinX2 must be non-negative")
    if args.model == "GP" and args.spinX2 % 2 != 0:
        raise ValueError("GP resonance spin must be integer")
    if args.model == "GP" and (not args.c_symmetry or not args.p_symmetry):
        raise ValueError("GP requires CP = [true,true], --no-c and --no-p are not supported")
    if mmax < 0:
        raise ValueError("--mmax must be non-negative")
    for label in (
        "mass",
        "width",
        "channel_mag",
        "channel_phase",
        "lambda_scale",
        "tol",
    ):
        value = float(getattr(args, label))
        if not math.isfinite(value):
            option = "Lambda" if label == "lambda_scale" else label.replace("_", "-")
            raise ValueError(f"--{option} must be finite")
    if args.mass < 0.0 or args.width < 0.0:
        raise ValueError("--mass and --width must be non-negative")
    if args.channel_mag < 0.0:
        raise ValueError("--channel-mag must be non-negative")
    if args.lambda_scale <= 0.0:
        raise ValueError("--Lambda must be positive")
    diphoton = tuple(args.fuse) == (22, 22)
    if args.model == "XP" and diphoton and args.spinX2 == 2:
        raise ValueError("XP gamma gamma fusion cannot produce a spin 1 resonance")
    if args.model == "XP" and diphoton and args.cparity != 1:
        raise ValueError("XP gamma gamma fusion requires C = +1")
    if diphoton and args.channel_mag != 1.0:
        raise ValueError("gamma gamma coupling magnitude is fixed by its partial width")
    if args.tol <= 0.0:
        raise ValueError("--tol must be positive")


# Parse command-line arguments
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Build covariant XP or analytic GP two exchange fusion couplings')
    parser.add_argument("--model", choices=("GP", "XP"), required=True)
    parser.add_argument("--fuse", nargs=2, type=int, metavar=("PDG1", "PDG2"), required=True)
    parser.add_argument(
        "--spinX2", type=int, required=True, help="produced resonance spin times two"
    )
    parser.add_argument("--P", dest="parity", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--C", dest="cparity", type=int, choices=(-1, 0, 1), required=True)
    parser.add_argument("--pdg", type=int, default=0, help="produced resonance PDG")
    parser.add_argument("--name", default="RES_generated", help="produced resonance name")
    parser.add_argument("--mass", type=float, default=0.0)
    parser.add_argument("--width", type=float, default=0.0)
    parser.add_argument(
        "--basis",
        choices=("g_ls", "helicity"),
        default=None,
        help="defaults to g_ls",
    )
    parser.add_argument(
        "--active-ls", nargs=2, metavar=("L", "S"), help="active LS row; default first allowed"
    )
    parser.add_argument(
        "--helicity-row",
        nargs=4,
        type=float,
        action="append",
        metavar=("M1", "M2", "MAG", "PHASE"),
        help="explicit GP helicity row, repeat for more rows",
    )
    parser.add_argument(
        "--channel-mag",
        type=float,
        default=1.0,
        help="absolute production coupling magnitude",
    )
    parser.add_argument(
        "--channel-phase", type=float, default=0.0, help="absolute production coupling phase"
    )
    parser.add_argument(
        "--Lambda",
        "--lambda-scale",
        dest="lambda_scale",
        type=float,
        default=1.0,
        help="LS momentum scale in GeV",
    )
    parser.add_argument("--tune-dir", default="modeldata/TUNE0")
    parser.add_argument("--mmax", type=int, default=None, help="override NUMERICS_REGGE.MMAX for GP")
    parser.add_argument(
        "--no-c", dest="c_symmetry", action="store_false", help="disable C-parity filtering"
    )
    parser.add_argument(
        "--no-p",
        dest="p_symmetry",
        action="store_false",
        help="disable parity/naturality filtering",
    )
    parser.add_argument(
        "--format",
        choices=("card", "json", "table"),
        default="card",
        help="compact card, expanded JSON, or derivation tables",
    )
    parser.add_argument("--tol", type=float, default=1e-12)
    return parser.parse_args()


# Execute the CLI
def main() -> int:
    args = parse_args()
    tune_dir = Path(args.tune_dir)
    particles = load_particles(tune_dir, args.model)
    default_mmax, regge = load_regge(tune_dir)
    mmax = default_mmax if args.mmax is None else args.mmax
    validate_args(args, mmax)
    fragment = build(args, particles, mmax, regge)
    if args.format == "table":
        print_table(args, fragment, particles, mmax)
    elif args.format == "json":
        print(json.dumps(fragment, indent=2, ensure_ascii=False))
    else:
        print(common.dumps(fragment))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"resonance_card_builder.py: {exc}", file=sys.stderr)
        raise SystemExit(1) from None
