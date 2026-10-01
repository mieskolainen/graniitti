#!/usr/bin/env python3
# Print symbolic Jacob-Wick LS and helicity basis conversions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
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

ORBITAL_LABELS = (
    "S",
    "P",
    "D",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "O",
    "Q",
    "R",
    "T",
    "U",
    "V",
    "W",
    "X",
    "Y",
    "Z",
)


@dataclass(frozen=True)
class LSRow:
    """LS basis row"""

    ell: int
    s: sp.Rational


@dataclass(frozen=True)
class HelicityRow:
    """Helicity basis row"""

    lambda1: sp.Rational
    lambda2: sp.Rational


# Parse an integer or rational spin label
def parse_spin(text: str) -> sp.Rational:
    value = sp.Rational(text)
    if 2 * value != int(2 * value) or value < 0:
        raise argparse.ArgumentTypeError(
            f"spin label should be a non-negative integer or half-integer: {text}"
        )
    return value


# Parse a parity quantum number
def parse_parity(text: str) -> int:
    value = int(text)
    if value not in {-1, 1}:
        raise argparse.ArgumentTypeError(f"parity should be +1 or -1: {text}")
    return value


# Compute compact labels for integer and half-integer values
def format_spin(value: sp.Rational) -> str:
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


# Compute a compact unsigned symbol token for spin labels
def spin_token(value: sp.Rational) -> str:
    if value.q == 1:
        return str(value.p)
    return f"{value.p}d{value.q}"


# Compute a compact signed integer token for magnetic labels
def int_token(value: int) -> str:
    if value == 0:
        return "0"
    if value < 0:
        return f"m{abs(value)}"
    return f"p{value}"


# Compute a standard helicity symbol token
def helicity_token(value: sp.Rational) -> str:
    if value == 0:
        return "0"
    if abs(value) in {sp.Rational(1, 2), sp.Rational(1, 1)}:
        return "-" if value < 0 else "+"
    sign = "-" if value < 0 else "+"
    return f"{sign}{spin_token(abs(value))}"


# Compute helicity or spin projections from -s to +s
def spin_projections(spin: sp.Rational) -> list[sp.Rational]:
    two_spin = int(2 * spin)
    return [sp.Rational(two_m, 2) for two_m in range(-two_spin, two_spin + 1, 2)]


# Compute true when an LS row passes the optional parity filter
def parity_allowed(ls_row: LSRow, parity: int | None, p1: int | None, p2: int | None) -> bool:
    if parity is None:
        return True
    assert p1 is not None and p2 is not None
    orbital_parity = 1 if ls_row.ell % 2 == 0 else -1
    return p1 * p2 * orbital_parity == parity


# Build the default angular-momentum-allowed LS basis
def ls_basis(
    J: sp.Rational,
    s1: sp.Rational,
    s2: sp.Rational,
    parity: int | None,
    p1: int | None,
    p2: int | None,
) -> list[LSRow]:
    if sp.Rational(J - s1 - s2).q != 1:
        return []
    rows = []
    for S in spin_projections(s1 + s2):
        if abs(s1 - s2) > S or s1 + s2 < S:
            continue
        l_min = int(math.ceil(float(abs(J - S))))
        l_max = int(math.floor(float(J + S)))
        for L in range(l_min, l_max + 1):
            row = LSRow(L, S)
            if parity_allowed(row, parity, p1, p2):
                rows.append(row)
    return rows


# Parse an explicit comma-separated L:S LS basis
def parse_ls_rows(
    spec: str,
    parity: int | None,
    p1: int | None,
    p2: int | None,
) -> list[LSRow]:
    rows = []
    seen = set()
    for item in spec.split(","):
        fields = item.strip().split(":")
        if len(fields) != 2:
            raise ValueError("--ls rows should use L:S,L:S syntax")
        row = LSRow(int(fields[0]), parse_spin(fields[1]))
        if row.ell < 0:
            raise ValueError("LS orbital angular momentum L should be non-negative")
        if row in seen:
            raise ValueError(f"duplicate LS row {item}")
        seen.add(row)
        if parity_allowed(row, parity, p1, p2):
            rows.append(row)
    return rows


# Build helicity rows in the same negative-to-positive order as GRANIITTI
def helicity_rows(J: sp.Rational, s1: sp.Rational, s2: sp.Rational) -> list[HelicityRow]:
    rows = []
    for lambda1 in spin_projections(s1):
        for lambda2 in spin_projections(s2):
            if abs(lambda1 - lambda2) <= J:
                rows.append(HelicityRow(lambda1, lambda2))
    return rows


# Compute the GRANIITTI Jacob-Wick LS to helicity coefficient
def jacob_wick(
    J: sp.Rational,
    s1: sp.Rational,
    s2: sp.Rational,
    ls_row: LSRow,
    h_row: HelicityRow,
) -> sp.Expr:
    helicity = h_row.lambda1 - h_row.lambda2
    norm = sp.sqrt(sp.Rational(2 * ls_row.ell + 1, 1) / (2 * J + 1))
    orbit_spin = clebsch_gordan(ls_row.ell, ls_row.s, J, 0, helicity, helicity)
    spin_helicity = clebsch_gordan(s1, s2, ls_row.s, h_row.lambda1, -h_row.lambda2, helicity)
    return sp.simplify(norm * orbit_spin * spin_helicity)


# Build the LS to helicity transformation matrix
def transform(
    J: sp.Rational,
    s1: sp.Rational,
    s2: sp.Rational,
    ls_rows: list[LSRow],
    h_rows: list[HelicityRow],
) -> sp.Matrix:
    validate_basis(J, s1, s2, ls_rows, h_rows)
    return sp.Matrix(
        [
            [jacob_wick(J, s1, s2, ls_row, h_row) for ls_row in ls_rows]
            for h_row in h_rows
        ]
    )


# Compute a symbol for one LS row
def alpha_symbol(ls_row: LSRow) -> sp.Symbol:
    return sp.Symbol(f"a_{ls_row.ell}{spin_token(ls_row.s)}")


# Compute a symbol for one helicity row
def helicity_symbol(h_row: HelicityRow) -> sp.Symbol:
    return sp.Symbol(f"H_{helicity_token(h_row.lambda1)}{helicity_token(h_row.lambda2)}")


# Compute a symbol for one real spherical harmonic
def harmonic_symbol(j_value: int, m_value: int) -> sp.Symbol:
    return sp.Symbol(f"R_{j_value}{int_token(m_value)}")


# Compute a symbol for one CEP real partial-wave coefficient
def cep_symbol(j_value: int, m_value: int) -> sp.Symbol:
    return sp.Symbol(f"C_{j_value}{int_token(m_value)}")


# Compute the spectroscopic letter for one orbital angular momentum
def orbital_label(l_value: int) -> str:
    if l_value < len(ORBITAL_LABELS):
        return ORBITAL_LABELS[l_value]
    return f"L{l_value}"


# Compute the spectroscopic label ^{2S+1}L_J
def wave_label(ls_row: LSRow, J: sp.Rational) -> str:
    multiplicity = int(2 * ls_row.s + 1)
    return f"^{multiplicity}{orbital_label(ls_row.ell)}_{format_spin(J)}"


# Format one expression for the requested output style
def format_expr(expr: sp.Expr, output_format: str) -> str:
    if output_format == "latex":
        return sp.latex(expr)
    if output_format == "plain":
        return sp.sstr(expr)
    return sp.pretty(expr)


# Format one exact expression for use inside a single-line table cell
def table_expr(expr: sp.Expr, output_format: str) -> str:
    simplified = sp.simplify(expr)
    if output_format == "latex":
        return sp.latex(simplified)
    return sp.sstr(simplified)


# Print equation with readable labels and exact symbolic expressions
def print_equation(
    lhs: str, expr: sp.Expr, output_format: str, *, simplify_expr: bool = True
) -> None:
    if simplify_expr:
        expr = sp.simplify(expr)
    rendered = format_expr(expr, output_format)
    if "\n" in rendered:
        print(f"{lhs} =")
        for line in rendered.splitlines():
            print(f"  {line}")
    else:
        print(f"{lhs} = {rendered}")


# Print labels and non-zero matrix coefficients for inspection
def print_basis(
    matrix: sp.Matrix, ls_rows: list[LSRow], h_rows: list[HelicityRow], output_format: str
) -> None:
    print("Basis")
    ls_table = [
        [
            index,
            row.ell,
            format_spin(row.s),
            alpha_symbol(row),
            f"^{int(2 * row.s + 1)}{orbital_label(row.ell)}",
        ]
        for index, row in enumerate(ls_rows)
    ]
    print(common.table(["index", "L", "S", "symbol", "LS term"], ls_table, right_align={0, 1, 2}))

    helicity_table = [
        [index, format_spin(row.lambda1), format_spin(row.lambda2), helicity_symbol(row)]
        for index, row in enumerate(h_rows)
    ]
    common.print_table("\nHelicity basis", ["index", "lambda1", "lambda2", "symbol"], helicity_table, right_align={0, 1, 2})

    coefficient_rows = []
    for i, h_row in enumerate(h_rows):
        coefficient_rows.append(
            [f"({format_spin(h_row.lambda1)},{format_spin(h_row.lambda2)})"]
            + [table_expr(matrix[i, j], output_format) for j in range(len(ls_rows))]
        )
    headers = ["(lambda1,lambda2)"] + [f"(L={row.ell},S={format_spin(row.s)})" for row in ls_rows]
    common.print_table(
        "\nJacob-Wick matrix M[helicity, LS]", headers, coefficient_rows, right_align=set(range(1, len(headers)))
    )


# Print the LS to helicity conversion
def print_helicities(
    matrix: sp.Matrix,
    ls_rows: list[LSRow],
    h_rows: list[HelicityRow],
    output_format: str,
) -> None:
    alpha = sp.Matrix([alpha_symbol(row) for row in ls_rows])
    converted = matrix * alpha

    print("")
    print("LS to helicity")
    for i, row in enumerate(h_rows):
        lhs = f"  H[{format_spin(row.lambda1)},{format_spin(row.lambda2)}]"
        print_equation(lhs, converted[i], output_format)


# Print the helicity to LS conversion from the exact pseudoinverse
def print_ls(
    matrix: sp.Matrix,
    ls_rows: list[LSRow],
    h_rows: list[HelicityRow],
    output_format: str,
) -> None:
    helicity = sp.Matrix([helicity_symbol(row) for row in h_rows])
    inverse = matrix.pinv()
    converted = inverse * helicity

    print("")
    print("Helicity to LS")
    for i, row in enumerate(ls_rows):
        lhs = f"  {alpha_symbol(row)}"
        print_equation(lhs, converted[i], output_format)


# Print spectroscopic labels in terms of compact LS and helicity amplitudes
def print_waves(
    matrix: sp.Matrix,
    J: sp.Rational,
    ls_rows: list[LSRow],
    h_rows: list[HelicityRow],
    output_format: str,
) -> None:
    helicity = sp.Matrix([helicity_symbol(row) for row in h_rows])
    converted = matrix.pinv() * helicity

    print("")
    print("Spectroscopic LS amplitudes")
    for i, row in enumerate(ls_rows):
        lhs = f"  {wave_label(row, J)} := {alpha_symbol(row)}"
        print_equation(lhs, converted[i], output_format)


# Compute true when a J value passes the requested CEP wave filter
def cep_wave_is_active(j_value: int, waves: str) -> bool:
    if waves == "even":
        return j_value % 2 == 0
    if waves == "odd":
        return j_value % 2 == 1
    return True


# Compute the expanded real tesseral harmonic in the associated-Legendre convention
def real_harmonic(
    j_value: int, m_value: int, theta: sp.Symbol, phi: sp.Symbol
) -> sp.Expr:
    m_abs = abs(m_value)
    norm = sp.sqrt(
        sp.Rational(2 * j_value + 1, 1)
        * sp.factorial(j_value - m_abs)
        / (4 * sp.pi * sp.factorial(j_value + m_abs))
    )
    legendre = sp.assoc_legendre(j_value, m_abs, sp.cos(theta))
    legendre = legendre.xreplace({sp.Abs(sp.sin(theta)): sp.sin(theta)})

    if m_value == 0:
        expr = norm * legendre
    elif m_value > 0:
        expr = sp.sqrt(2) * norm * legendre * sp.cos(m_abs * phi)
    else:
        expr = sp.sqrt(2) * norm * legendre * sp.sin(m_abs * phi)
    expr = sp.simplify(expr)
    expr = expr.replace(lambda x: x == sp.Abs(sp.sin(theta)), lambda x: sp.sin(theta))
    return sp.simplify(expr)


# Print the generic CEP pion-pair real partial-wave expansion
def print_cep(jmax: int, waves: str, output_format: str) -> None:
    if jmax < 0:
        raise ValueError("--jmax should be non-negative")

    theta_pi, phi_pi = sp.symbols("theta_pi phi_pi", real=True)
    active_j = [j_value for j_value in range(jmax + 1) if cep_wave_is_active(j_value, waves)]
    if not active_j:
        raise ValueError("empty CEP wave set after --waves filter")

    amplitude = 0
    for j_value in active_j:
        for m_value in range(-j_value, j_value + 1):
            amplitude += cep_symbol(j_value, m_value) * harmonic_symbol(
                j_value, m_value
            )

    print("Central pp -> pp pi pi real partial-wave expansion")
    print("Process: p(lambda1) p(lambda2) -> p(lambda3) p(lambda4) pi pi")
    print("Angles: theta_pi, phi_pi are evaluated in the pi-pi rest frame")
    print(
        'C_JM depends on proton helicities, t1, t2, Delta_phi_pp and m_pipi'
    )
    print("")
    print("Measured-event amplitude")
    print_equation("  A_{lambda3 lambda4; lambda1 lambda2}", amplitude, output_format)

    print("")
    print("Coefficient factorization")
    print(
        "  C_JM := BW_J(m_pipi) g_{J->pipi} V_JM(lambda1,lambda2,lambda3,lambda4;t1,t2,Delta_phi_pp)"
    )
    print("  V_JM := sum_{m1,m2} B_up(lambda1->lambda3,m1;t1) B_dn(lambda2->lambda4,m2;t2)")
    print("          * exp(i m1 Delta_phi_pp) delta_{M,m1-m2} P_JM(m1,m2)")
    print("  P_JM := sum_{L,S} a_prod[J,L,S] sqrt((2L+1)/(2J+1))")
    print("          * <L 0, S M | J M> <j1 m1, j2 -m2 | S M>")
    print(
        '  X_J -> pi pi: L=J, S=0, decay basis R_JM(theta_pi,phi_pi)'
    )

    print("")
    print("Real spherical harmonics")
    print('  Condon-Shortley convention: m>0 cosine, m<0 sine')
    harmonic_rows = []
    for j_value in active_j:
        for m_value in range(-j_value, j_value + 1):
            expression = real_harmonic(j_value, m_value, theta_pi, phi_pi)
            harmonic_rows.append(
                [
                    j_value,
                    m_value,
                    harmonic_symbol(j_value, m_value),
                    table_expr(expression, output_format),
                ]
            )
    print(
        common.table(
            ["J", "M", "symbol", "R_JM(theta_pi,phi_pi)"],
            harmonic_rows,
            right_align={0, 1},
            group_by=0,
        )
    )


# Validate the angular momentum triangles and spin projections of the requested basis
def validate_basis(
    J: sp.Rational,
    s1: sp.Rational,
    s2: sp.Rational,
    ls_rows: list[LSRow],
    h_rows: list[HelicityRow],
) -> None:
    if sp.Rational(J - s1 - s2).q != 1:
        raise ValueError("mother and daughter spins must differ by an integer")
    if not ls_rows:
        raise ValueError("empty LS basis after filters")
    if not h_rows:
        raise ValueError("empty helicity basis")
    allowed_ls = set(ls_basis(J, s1, s2, None, None, None))
    for row in ls_rows:
        if row not in allowed_ls:
            raise ValueError(f"forbidden LS row L={row.ell}, S={row.s} for J={J}, s1={s1}, s2={s2}")
    allowed_helicity = set(helicity_rows(J, s1, s2))
    if any(row not in allowed_helicity for row in h_rows):
        raise ValueError("forbidden helicity row for the requested spins")
    if len(set(ls_rows)) != len(ls_rows) or len(set(h_rows)) != len(h_rows):
        raise ValueError("duplicate rows in the requested basis")


# Build the command line parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Compute exact Jacob-Wick LS and helicity conversions')
    parser.add_argument(
        "--mode", choices=("ls", "cep-pipi"), default="ls", help="Expression family to print"
    )
    parser.add_argument("--J", type=parse_spin, help="Mother spin, e.g. 0, 1, 3/2")
    parser.add_argument("--s1", type=parse_spin, help="First daughter spin")
    parser.add_argument("--s2", type=parse_spin, help="Second daughter spin")
    parser.add_argument("--ls", help="Optional comma-separated LS rows as L:S,L:S")
    parser.add_argument(
        "--parity", type=parse_parity, help="Optional mother parity filter, +1 or -1"
    )
    parser.add_argument("--p1", type=parse_parity, help="First daughter parity for --parity")
    parser.add_argument("--p2", type=parse_parity, help="Second daughter parity for --parity")
    parser.add_argument(
        "--jmax", type=int, default=2, help="Maximum central pi-pi spin J for --mode cep-pipi"
    )
    parser.add_argument(
        "--waves",
        choices=("all", "even", "odd"),
        default="all",
        help="J-wave filter for --mode cep-pipi",
    )
    parser.add_argument(
        "--format", choices=("pretty", "plain", "latex"), default="pretty", help="Output format"
    )
    return parser


# Run the symbolic conversion printer
def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.mode == "cep-pipi":
        print_cep(args.jmax, args.waves, args.format)
        return 0

    if args.J is None or args.s1 is None or args.s2 is None:
        parser.error("--mode ls requires --J, --s1 and --s2")

    parity_options = (args.parity, args.p1, args.p2)
    if any(value is not None for value in parity_options) and any(
        value is None for value in parity_options
    ):
        parser.error("--parity requires both --p1 and --p2, and --p1/--p2 require --parity")

    if args.ls:
        ls_rows = parse_ls_rows(args.ls, args.parity, args.p1, args.p2)
    else:
        ls_rows = ls_basis(args.J, args.s1, args.s2, args.parity, args.p1, args.p2)
    h_rows = helicity_rows(args.J, args.s1, args.s2)

    matrix = transform(args.J, args.s1, args.s2, ls_rows, h_rows)

    print("GRANIITTI Jacob-Wick LS <-> helicity conversion")
    print(f"J={format_spin(args.J)}, s1={format_spin(args.s1)}, s2={format_spin(args.s2)}")
    print("Convention: lambda = lambda1 - lambda2")
    print(
        "Coefficient: sqrt((2L+1)/(2J+1)) * <L 0, S lambda | J lambda> * <s1 lambda1, s2 -lambda2 | S lambda>"
    )
    print("")
    print_basis(matrix, ls_rows, h_rows, args.format)
    print_helicities(matrix, ls_rows, h_rows, args.format)
    print_ls(matrix, ls_rows, h_rows, args.format)
    print_waves(matrix, args.J, ls_rows, h_rows, args.format)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (TypeError, ValueError) as exc:
        print(f"ls_helicity_conversion.py: {exc}", file=sys.stderr)
        raise SystemExit(1) from None
