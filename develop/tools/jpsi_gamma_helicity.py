#!/usr/bin/env python3
# Convert chi_cJ -> J/psi gamma multipoles to GRANIITTI helicity rows
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
#
# [REFERENCE: Artuso et al., arXiv:0910.0046, Eqs. (4)-(6)]
# [REFERENCE: BESIII, arXiv:1701.01197]
# [REFERENCE: Particle Data Group 2026 chi_c1 and chi_c2 listings]

import argparse
import json
import math
from dataclasses import dataclass

PDG_2026_MULTIPOLES = {
    0: (0.0, 0.0),
    1: (-0.067, 0.0),
    2: (-0.110, -0.003),
}


@dataclass(frozen=True)
class Conversion:
    """Normalized multipoles and helicity rows"""

    spin: int
    multipoles: tuple[float, ...]
    helicities: tuple[float, ...]
    rows: tuple[tuple[int, int, float, float], ...]


# Compute normalized E1, M2 and E3 amplitudes with the positive-E1 convention
def normalize(spin: int, m2: float = 0.0, e3: float = 0.0) -> tuple[float, ...]:
    if spin not in (0, 1, 2):
        raise ValueError("spin must be 0, 1 or 2")
    if not math.isfinite(m2) or not math.isfinite(e3):
        raise ValueError("multipole amplitudes must be finite")
    if spin == 0:
        if abs(m2) > 0.0 or abs(e3) > 0.0:
            raise ValueError("chi_c0 -> J/psi gamma has only the E1 multipole")
        return (1.0,)
    if spin == 1 and abs(e3) > 0.0:
        raise ValueError("chi_c1 -> J/psi gamma does not allow an E3 multipole")

    norm2 = m2 * m2 + e3 * e3
    if norm2 >= 1.0:
        raise ValueError("normalized M2 and E3 amplitudes require M2^2 + E3^2 < 1")
    e1 = math.sqrt(1.0 - norm2)
    return (e1, m2) if spin == 1 else (e1, m2, e3)


# Transform normalized multipoles into the standard real A_nu helicity amplitudes
def to_helicity(spin: int, multipoles: tuple[float, ...]) -> tuple[float, ...]:
    if spin == 0:
        if len(multipoles) != 1:
            raise ValueError("spin zero expects one E1 amplitude")
        return multipoles
    if spin == 1:
        if len(multipoles) != 2:
            raise ValueError("spin one expects E1 and M2 amplitudes")
        e1, m2 = multipoles
        return ((e1 + m2) / math.sqrt(2.0), (e1 - m2) / math.sqrt(2.0))
    if spin == 2:
        if len(multipoles) != 3:
            raise ValueError("spin two expects E1, M2 and E3 amplitudes")
        e1, m2, e3 = multipoles
        return (
            math.sqrt(1.0 / 10.0) * e1 + math.sqrt(1.0 / 2.0) * m2 + math.sqrt(2.0 / 5.0) * e3,
            math.sqrt(3.0 / 10.0) * e1 + math.sqrt(1.0 / 6.0) * m2 - math.sqrt(8.0 / 15.0) * e3,
            math.sqrt(3.0 / 5.0) * e1 - math.sqrt(1.0 / 3.0) * m2 + math.sqrt(1.0 / 15.0) * e3,
        )
    raise ValueError("spin must be 0, 1 or 2")


# Compute compact GRANIITTI rows normalized to one reference helicity amplitude
def card_rows(
    spin: int, helicities: tuple[float, ...]
) -> tuple[tuple[int, int, float, float], ...]:
    row_map = {
        0: (((-1, -1), 0),),
        1: (((-1, -1), 0), ((0, -1), 1)),
        2: (((-1, -1), 0), ((-1, 1), 2), ((0, -1), 1)),
    }[spin]
    reference_index = 2 if spin == 2 else 0
    reference = helicities[reference_index]
    if abs(reference) < 1.0e-15:
        raise ValueError("the selected reference helicity amplitude is zero")

    rows = []
    for (lambda1, lambda2), amplitude_index in row_map:
        amplitude = helicities[amplitude_index]
        ratio = amplitude / reference
        phase = -math.pi if ratio < 0.0 else 0.0
        rows.append((lambda1, lambda2, abs(ratio), phase))
    return tuple(rows)


# Convert one chi_c spin and normalized higher multipoles
def convert(spin: int, m2: float = 0.0, e3: float = 0.0) -> Conversion:
    multipoles = normalize(spin, m2, e3)
    helicities = to_helicity(spin, multipoles)
    norm2 = sum(value * value for value in helicities)
    if not math.isclose(norm2, 1.0, rel_tol=0.0, abs_tol=2.0e-15):
        raise RuntimeError(f"multipole transformation is not orthogonal: norm^2={norm2}")
    return Conversion(
        spin=spin,
        multipoles=multipoles,
        helicities=helicities,
        rows=card_rows(spin, helicities),
    )


# Compute conversion as a JSON-serializable dictionary
def as_dict(result: Conversion) -> dict:
    names = ("E1", "M2", "E3")[: len(result.multipoles)]
    return {
        "spin": result.spin,
        "multipoles": dict(zip(names, result.multipoles, strict=False)),
        "helicities": {f"A{index}": value for index, value in enumerate(result.helicities)},
        "helicity": [list(row) for row in result.rows],
    }


# Print conversion in a compact human-readable form
def print_table(result: Conversion) -> None:
    payload = as_dict(result)
    multipoles = ", ".join(f"{name}={value:.16g}" for name, value in payload["multipoles"].items())
    helicities = ", ".join(f"{name}={value:.16g}" for name, value in payload["helicities"].items())
    print(f"chi_c{result.spin} -> J/psi gamma")
    print(f"  multipoles: {multipoles}")
    print(f"  helicities: {helicities}")
    print('  "helicity" : ' + json.dumps(payload["helicity"], separators=(",", ":")))


# Build the command-line parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Convert chi_cJ -> J/psi gamma multipoles to helicity rows')
    parser.add_argument("--spin", type=int, choices=(0, 1, 2))
    parser.add_argument(
        "--preset",
        choices=("pdg2026", "pure-e1"),
        default="pdg2026",
        help="input multipoles used when --m2 or --e3 is not supplied",
    )
    parser.add_argument("--m2", type=float, help="normalized M2 amplitude")
    parser.add_argument("--e3", type=float, help="normalized E3 amplitude")
    parser.add_argument("--json", action="store_true", help="print machine-readable JSON")
    return parser


# Resolve command-line inputs and print one or all chi_cJ conversions
def main() -> None:
    args = build_parser().parse_args()
    spins = (args.spin,) if args.spin is not None else (0, 1, 2)
    results = []
    for spin in spins:
        preset_m2, preset_e3 = PDG_2026_MULTIPOLES[spin] if args.preset == "pdg2026" else (0.0, 0.0)
        m2 = preset_m2 if args.m2 is None else args.m2
        e3 = preset_e3 if args.e3 is None else args.e3
        results.append(convert(spin, m2, e3))

    if args.json:
        print(json.dumps([as_dict(result) for result in results], indent=2))
        return
    for index, result in enumerate(results):
        if index:
            print()
        print_table(result)


if __name__ == "__main__":
    main()
