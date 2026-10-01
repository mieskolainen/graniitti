#!/usr/bin/env python3
# Derive pure central |Jz| sectors at equal integer exchange poles
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import sympy as sp
from sympy.physics.wigner import clebsch_gordan

# Resolve the development package when this file is executed directly
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from develop.tools.lib import common


@dataclass(frozen=True)
class Model:
    """Central fusion quantum numbers"""

    spin: int
    parity: int
    mmax: int
    exchange_naturality: int
    identical: bool
    pole_spin: int


@dataclass(frozen=True)
class LSRow:
    """LS quantum numbers"""

    ell: int
    spin: int


Pair = tuple[int, int]
Vector = dict[Pair, sp.Expr]


# Compute a validated sign quantum number
def parse_sign(text: str) -> int:
    value = int(text)
    if value not in {-1, 1}:
        raise argparse.ArgumentTypeError("sign should be +1 or -1")
    return value


# Compute a comma-separated integer list
def parse_spins(text: str) -> list[int]:
    values = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        value = int(item)
        if value < 0:
            raise argparse.ArgumentTypeError("spin values should be non-negative integers")
        values.append(value)
    if not values:
        raise argparse.ArgumentTypeError("at least one spin value is required")
    return values


# Compute (-1)^n for integer n
def phase_sign(n: int) -> int:
    return 1 if n % 2 == 0 else -1


# Compute the Clebsch-Gordan coefficient in the GRANIITTI argument order
def cg(j1: int, j2: int, m1: int, m2: int, j: int, m: int) -> sp.Expr:
    return sp.simplify(clebsch_gordan(j1, j2, j, m1, m2, m))


# Compute the analytic central exchange-helicity pairs kept by MMAX
def helicity_pairs(model: Model) -> list[Pair]:
    out = []
    limit = min(model.mmax, model.pole_spin)
    for m1 in range(-limit, limit + 1):
        for m2 in range(-limit, limit + 1):
            if abs(m1 - m2) <= model.spin:
                out.append((m1, m2))
    return out


# Compute true if an LS row passes GP PP selection rules
def ls_allowed(model: Model, ell: int, spin: int) -> bool:
    if spin < 0 or ell < 0:
        return False
    if abs(ell - spin) > model.spin or model.spin > ell + spin:
        return False
    if model.exchange_naturality * phase_sign(ell) != model.parity:
        return False
    return not (model.identical and (ell - spin) % 2 != 0)


# Compute all LS rows allowed by the fixed exchange spins
def allowed_ls_rows(model: Model) -> list[LSRow]:
    rows = []
    max_spin = 2 * model.pole_spin
    for spin in range(0, max_spin + 1):
        for ell in range(0, model.spin + spin + 1):
            if ls_allowed(model, ell, spin):
                rows.append(LSRow(ell, spin))
    return rows


# Compute pole LS image component up to a nonzero column normalization
def ls_component(model: Model, row: LSRow, pair: Pair) -> sp.Expr:
    m1, m2 = pair
    mu = m1 - m2
    if abs(mu) > model.spin or abs(mu) > row.spin:
        return sp.Integer(0)
    norm = sp.sqrt(sp.Rational(2 * row.ell + 1, 2 * model.spin + 1))
    orbital_spin = cg(row.ell, row.spin, 0, mu, model.spin, mu)
    helicity_spin = cg(model.pole_spin, model.pole_spin, m1, -m2, row.spin, mu)
    return sp.simplify(norm * orbital_spin * helicity_spin)


# Compute the exact LS image matrix and raw column vectors
def ls_image(
    model: Model, rows: list[LSRow], pairs: list[Pair]
) -> tuple[sp.Matrix, list[list[sp.Expr]]]:
    columns = []
    for row in rows:
        columns.append([ls_component(model, row, pair) for pair in pairs])
    matrix = (
        sp.Matrix.hstack(*(sp.Matrix(col) for col in columns))
        if columns
        else sp.zeros(len(pairs), 0)
    )
    return matrix, columns


# Compute the exact squared norm of a sparse vector
def norm_sq(vector: Vector) -> sp.Expr:
    return sp.simplify(sum(sp.conjugate(value) * value for value in vector.values()))


# Compute a normalized sparse vector with a positive real seed convention
def normalize(vector: Vector) -> Vector:
    norm2 = norm_sq(vector)
    if norm2 == 0:
        raise ValueError("cannot normalize a zero vector")
    norm = sp.sqrt(norm2)
    return {pair: sp.simplify(value / norm) for pair, value in vector.items()}


# Compute true if a sparse vector lies in the LS image span
def in_span(matrix: sp.Matrix, pairs: list[Pair], vector: Vector) -> bool:
    dense = sp.Matrix([vector.get(pair, sp.Integer(0)) for pair in pairs])
    return matrix.row_join(dense).rank() == matrix.rank()


# Compute the symmetry phase between two H coordinates from the LS image
def symmetry_phase(
    columns: list[list[sp.Expr]],
    pair_index: dict[Pair, int],
    source: Pair,
    target: Pair,
) -> sp.Expr | None:
    i = pair_index[source]
    j = pair_index[target]
    phase = None
    for column in columns:
        a = sp.simplify(column[i])
        b = sp.simplify(column[j])
        if a == 0 and b == 0:
            continue
        if a == 0 or b == 0:
            return None
        ratio = sp.simplify(b / a)
        if ratio == 0:
            return None
        unit = sp.simplify(ratio / sp.sqrt(sp.conjugate(ratio) * ratio))
        if phase is None:
            phase = unit
        elif sp.simplify(phase - unit) != 0:
            return None
    return phase


# Insert a sparse-vector entry and reject incompatible orbit constraints
def insert_entry(vector: Vector, pair: Pair, value: sp.Expr) -> bool:
    value = sp.simplify(value)
    if value == 0:
        return False
    if pair in vector:
        if sp.simplify(vector[pair] - value) != 0:
            raise ValueError(f"inconsistent orbit constraint for pair {pair}")
        return False
    vector[pair] = value
    return True


# Expand one seed by the LS-derived parity and identical-exchange phases
def expand_orbit(
    seed: Pair,
    columns: list[list[sp.Expr]],
    pair_index: dict[Pair, int],
    identical: bool,
) -> Vector:
    vector: Vector = {}
    queue: list[Pair] = []
    insert_entry(vector, seed, sp.Integer(1))
    queue.append(seed)

    while queue:
        pair = queue.pop(0)
        value = vector[pair]

        parity_pair = (-pair[0], -pair[1])
        if parity_pair in pair_index:
            phase = symmetry_phase(columns, pair_index, pair, parity_pair)
            if phase is not None and insert_entry(vector, parity_pair, phase * value):
                queue.append(parity_pair)

        exchange_pair = (pair[1], pair[0])
        if identical and exchange_pair in pair_index:
            phase = symmetry_phase(columns, pair_index, pair, exchange_pair)
            if phase is not None and insert_entry(vector, exchange_pair, phase * value):
                queue.append(exchange_pair)

    return normalize(vector)


# Compute candidate seed pairs ordered by compactness
def candidate_seeds(model: Model, sector: int) -> list[Pair]:
    seeds = []
    for pair in helicity_pairs(model):
        mu = pair[0] - pair[1]
        if (sector == 0 and mu == 0) or (sector > 0 and abs(mu) == sector):
            seeds.append(pair)
    return sorted(
        seeds, key=lambda p: (max(abs(p[0]), abs(p[1])), abs(p[0]) + abs(p[1]), p[0], p[1])
    )


# Compute sparse pure-|Jz| representative if it exists
def sector_vector(
    model: Model,
    sector: int,
    matrix: sp.Matrix,
    columns: list[list[sp.Expr]],
    pairs: list[Pair],
) -> Vector | None:
    pair_index = {pair: i for i, pair in enumerate(pairs)}
    for seed in candidate_seeds(model, sector):
        try:
            vector = expand_orbit(seed, columns, pair_index, model.identical)
        except ValueError:
            continue
        if any(abs(pair[0] - pair[1]) != sector for pair in vector):
            continue
        if in_span(matrix, pairs, vector):
            return vector
    return None


# Compute magnitude and phase for one complex table entry
def mag_phase(value: sp.Expr) -> tuple[float, float]:
    z = complex(sp.N(value, 18))
    if abs(z.real) < 1e-15:
        z = complex(0.0, z.imag)
    if abs(z.imag) < 1e-15:
        z = complex(z.real, 0.0)
    mag = abs(z)
    phase = common.phase(math.atan2(z.imag, z.real))
    if abs(phase) < 1e-15:
        phase = 0.0
    return mag, phase


# Compute a compact matrix string with rows and columns ordered by m
def format_matrix(vector: Vector, mmax: int) -> str:
    rows = []
    for m1 in range(-mmax, mmax + 1):
        fields = []
        for m2 in range(-mmax, mmax + 1):
            value = vector.get((m1, m2), sp.Integer(0))
            z = complex(sp.N(value, 12))
            if abs(z.imag) < 1e-12:
                fields.append(f"{z.real: .6g}")
            else:
                fields.append(f"{z.real: .6g}{z.imag:+.6g}i")
        rows.append("[" + ", ".join(fields) + "]")
    return "\n".join(rows)


# Build all pure-|Jz| representatives for one model
def derive(model: Model) -> dict[str, object]:
    rows = allowed_ls_rows(model)
    pairs = helicity_pairs(model)
    matrix, columns = ls_image(model, rows, pairs)
    jp = f"{model.spin}{'+' if model.parity > 0 else '-'}"
    sectors = []
    for sector in range(0, model.spin + 1):
        vector = sector_vector(model, sector, matrix, columns, pairs)
        entries = []
        if vector is not None:
            for pair, value in sorted(vector.items()):
                magnitude, phase = mag_phase(value)
                entries.append([pair[0], pair[1], magnitude, phase])
        sectors.append(
            {
                "abs_Jz": sector,
                "found": vector is not None,
                "helicity": entries,
                "matrix": format_matrix(vector, model.mmax) if vector is not None else None,
            }
        )
    return {
        "J": model.spin,
        "P": model.parity,
        "J_P": jp,
        "MMAX": model.mmax,
        "pole_spin": model.pole_spin,
        "exchange_naturality_product": model.exchange_naturality,
        "identical_exchange": model.identical,
        "allowed_ls": [[row.ell, row.spin] for row in rows],
        "ls_image_rank": matrix.rank(),
        "sectors": sectors,
    }


# Print model summaries and sparse helicity representatives as aligned tables
def print_tables(results: list[dict[str, object]], show_matrices: bool) -> None:
    summary_rows = []
    for result in results:
        ls_rows = ", ".join(f"({row[0]},{row[1]})" for row in result["allowed_ls"])
        found = sum(1 for sector in result["sectors"] if sector["found"])
        summary_rows.append(
            [
                result["J_P"],
                result["MMAX"],
                result["exchange_naturality_product"],
                "yes" if result["identical_exchange"] else "no",
                ls_rows or "none",
                result["ls_image_rank"],
                f"{found}/{len(result['sectors'])}",
            ]
        )
    common.print_table(
        "GP pure-|Jz| helicity representatives",
        ["J^P", "MMAX", "eta1*eta2", "identical", "allowed (L,S)", "LS rank", "sectors"],
        summary_rows,
        right_align={1, 2, 5, 6},
    )

    helicity_rows = []
    for result in results:
        for sector in result["sectors"]:
            label = "0" if sector["abs_Jz"] == 0 else f"+/-{sector['abs_Jz']}"
            if not sector["found"]:
                helicity_rows.append([result["J_P"], label, "-", "-", "not found", "-"])
                continue
            for m1, m2, magnitude, phase in sector["helicity"]:
                helicity_rows.append(
                    [result["J_P"], label, m1, m2, f"{magnitude:.12g}", f"{phase:.12g}"]
                )
    common.print_table(
        "\nSparse normalized helicity rows",
        ["J^P", "|Jz|", "m1", "m2", "magnitude", "phase [rad]"],
        helicity_rows,
        right_align={1, 2, 3, 4, 5},
        group_by=0,
    )

    if show_matrices:
        for result in results:
            for sector in result["sectors"]:
                if not sector["found"]:
                    continue
                label = "0" if sector["abs_Jz"] == 0 else f"+/-{sector['abs_Jz']}"
                print(f"\nJ^P={result['J_P']}, |Jz|={label}; matrix rows/columns m=-MMAX...+MMAX")
                print(sector["matrix"])


# Compute independent GP helicity couplings with symmetry partners left implicit
def card_rows(rows: list[list[int | float]], identical: bool) -> list[list[int | float]]:
    independent = []
    covered = set()
    for row in rows:
        m1, m2 = int(row[0]), int(row[1])
        orbit = [(m1, m2), (-m1, -m2)]
        if identical:
            orbit.extend([(m2, m1), (-m2, -m1)])
        key = min(orbit)
        if key not in covered:
            independent.append(list(row))
            covered.add(key)
    return independent


# Compute compact GP card alternatives with independent symmetry orbits
def cards(results: list[dict[str, object]]) -> dict[str, object]:
    cards = {}
    for result in results:
        cards[result["J_P"]] = {
            f"abs_Jz_{sector['abs_Jz']}": {
                "basis": "helicity",
                "CP": [True, True],
                "helicity": card_rows(sector["helicity"], result["identical_exchange"]),
            }
            for sector in result["sectors"]
            if sector["found"]
        }
    return cards


# Build the command-line parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Derive pure |Jz| helicity sectors at exchange poles')
    parser.add_argument(
        "--J", default="0,1,2", type=parse_spins, help="comma-separated integer central spins"
    )
    parser.add_argument("--parity", default=1, type=parse_sign, help="central parity")
    parser.add_argument("--pole-spin", required=True, type=int, help="common integer exchange pole spin")
    parser.add_argument(
        "--mmax",
        "--MMAX",
        dest="mmax",
        default=2,
        type=int,
        help="analytic exchange-helicity cutoff",
    )
    parser.add_argument(
        "--exchange-naturality",
        default=1,
        type=parse_sign,
        help="product eta1*eta2, +1 for Pomeron-Pomeron",
    )
    parser.add_argument(
        "--non-identical", action="store_true", help="disable identical-exchange symmetry"
    )
    parser.add_argument("--format", choices=("table", "json", "cards"), default="table")
    parser.add_argument(
        "--matrix", action="store_true", help="include expanded matrices in table output"
    )
    return parser


# Program entry point
def main() -> int:
    args = build_parser().parse_args()
    if args.mmax < 0:
        raise ValueError("--mmax should be non-negative")
    if args.pole_spin < 0:
        raise ValueError("--pole-spin should be non-negative")
    results = []
    for spin in args.J:
        model = Model(
            spin=spin,
            parity=args.parity,
            mmax=args.mmax,
            exchange_naturality=args.exchange_naturality,
            identical=not args.non_identical,
            pole_spin=args.pole_spin,
        )
        results.append(derive(model))
    if args.format == "table":
        print_tables(results, args.matrix)
    elif args.format == "json":
        print(json.dumps({"models": results}, indent=2, ensure_ascii=False))
    else:
        print(common.dumps(cards(results)))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (TypeError, ValueError) as exc:
        print(f"analytic_Jz.py: {exc}", file=sys.stderr)
        raise SystemExit(1) from None
