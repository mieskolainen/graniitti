# Shared particle metadata and pole spin algebra
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from typing import Any

from . import common


@dataclass(frozen=True)
class Particle:
    pdg: int
    name: str
    spin2: int
    parity: int
    cparity: int
    mass: float = 0.0
    width: float = 0.0
    charge3: int = 0
    isospin2: int = 0
    gparity: int = 0


@dataclass(frozen=True)
class ReggeTable:
    groups: list[list[int]]
    tau: list[int]
    pole_spin: list[int]


# Compute true if x is integral within numerical tolerance
def is_integer(x: float, tol: float = 1e-9) -> bool:
    return abs(x - round(x)) < tol


# Compute a JSON-friendly integer when a spin label is integral
def clean_number(x: float) -> int | float:
    if is_integer(x):
        return int(round(x))
    return round(x, 15)


# Compute n! generalized to integer and half-integer arguments
def fact(x: float) -> float:
    if x < -1e-9:
        return 0.0
    return math.gamma(x + 1.0)


# Evaluate the Condon-Shortley Clebsch-Gordan coefficient
def clebsch_gordan(j1: float, j2: float, m1: float, m2: float, j: float, m: float) -> float:
    if any(not math.isfinite(x) for x in (j1, j2, j, m1, m2, m)):
        raise ValueError("Clebsch-Gordan quantum numbers must be finite")
    if (any(x < 0 or not is_integer(2 * x) for x in (j1, j2, j))
            or any(not is_integer(x) for x in (j1 + m1, j2 + m2, j + m, j1 + j2 + j))
            or abs(m1 + m2 - m) > 1e-9):
        return 0.0
    if j < abs(j1 - j2) - 1e-9 or j > j1 + j2 + 1e-9:
        return 0.0
    for value in (j1 + m1, j1 - m1, j2 + m2, j2 - m2, j + m, j - m):
        if value < -1e-9:
            return 0.0

    pref = math.sqrt(
        (2.0 * j + 1.0)
        * fact(j1 + j2 - j)
        * fact(j1 - j2 + j)
        * fact(-j1 + j2 + j)
        / fact(j1 + j2 + j + 1.0)
    )
    pref *= math.sqrt(
        fact(j1 + m1) * fact(j1 - m1) * fact(j2 + m2) * fact(j2 - m2) * fact(j + m) * fact(j - m)
    )

    k_min = math.ceil(max(0.0, j2 - j - m1, j1 - j + m2) - 1e-9)
    k_max = math.floor(min(j1 + j2 - j, j1 - m1, j2 + m2) + 1e-9)
    total = 0.0
    for k in range(int(k_min), int(k_max) + 1):
        denom = (
            fact(k)
            * fact(j1 + j2 - j - k)
            * fact(j1 - m1 - k)
            * fact(j2 + m2 - k)
            * fact(j - j2 + m1 + k)
            * fact(j - j1 - m2 + k)
        )
        if denom != 0.0:
            total += ((-1.0) ** k) / denom
    return pref * total


# Compute magnitude and phase with small numerical noise removed
def mag_phase(re: float, im: float, tol: float) -> tuple[float, float]:
    mag = math.hypot(re, im)
    if mag <= tol:
        return 0.0, 0.0
    phase = common.phase(math.atan2(im, re))
    if abs(phase) < 1e-12:
        phase = 0.0
    return mag, phase


# Compute a particle table with common built-ins used in steering cards
def particles() -> dict[int, Particle]:
    rows = {
        22: Particle(22, "gamma", 2, -1, -1),
        111: Particle(111, "pi0", 0, -1, 1),
        211: Particle(211, "pi+", 0, -1, 0),
        311: Particle(311, "K0", 0, -1, 0),
        321: Particle(321, "K+", 0, -1, 0),
        2112: Particle(2112, "n0", 1, 1, 0),
        2212: Particle(2212, "p+", 1, 1, 0),
        113: Particle(113, "rho0", 2, -1, -1),
        221: Particle(221, "eta", 0, -1, 1),
        225: Particle(225, "f2_1270", 4, 1, 1),
        331: Particle(331, "eta_prime", 0, -1, 1),
        333: Particle(333, "phi", 2, -1, -1),
        335: Particle(335, "f2_1525", 4, 1, 1),
        9000221: Particle(9000221, "f0", 0, 1, 1),
    }
    out = dict(rows)
    for pdg in (211, 311, 321):
        row = rows[pdg]
        out[-pdg] = Particle(-pdg, row.name.replace("+", "-"), row.spin2, row.parity, row.cparity)
    for pdg in (2112, 2212):
        row = rows[pdg]
        out[-pdg] = Particle(-pdg, row.name.replace("+", "-"), row.spin2, -row.parity, row.cparity)
    return out


# Compute true when a particle has a defined C-parity eigenvalue
def has_cparity(particle: Particle) -> bool:
    return particle.cparity in (-1, 1)


# Compute the C parity of a two-body state when it is defined
def cparity(p1: Particle, p2: Particle, l_value: int, s_value: float) -> int | None:
    if p1.pdg == -p2.pdg:
        exponent = l_value + s_value
        if not is_integer(exponent):
            return None
        return 1 if int(round(exponent)) % 2 == 0 else -1
    if has_cparity(p1) and has_cparity(p2):
        return p1.cparity * p2.cparity
    return None


# Compute true if the two-body C-parity rule accepts one LS row
def c_allowed(
    mother: Particle, p1: Particle, p2: Particle, l_value: int, s_value: float, c_symmetry: bool
) -> bool:
    if not c_symmetry or not has_cparity(mother):
        return True
    c_total = cparity(p1, p2, l_value, s_value)
    return c_total is None or c_total == mother.cparity


# Compute whether one LS row has nonzero physical pole helicity support
def ls_supported(
    mother: Particle,
    p1: Particle,
    p2: Particle,
    l_value: int,
    s_value: float,
) -> bool:
    if mother.spin2 % 2 != 0 or any(particle.spin2 < 0 for particle in (mother, p1, p2)):
        return False
    if (p1.pdg == 22 and p1.spin2 != 2) or (p2.pdg == 22 and p2.spin2 != 2):
        return False
    if not is_integer(s_value):
        return False

    j_value = 0.5 * mother.spin2
    j1 = p1.spin2 / 2.0
    j2 = p2.spin2 / 2.0
    helicities1 = (-1, 1) if p1.pdg == 22 else tuple(-j1 + i for i in range(p1.spin2 + 1))
    helicities2 = (-1, 1) if p2.pdg == 22 else tuple(-j2 + i for i in range(p2.spin2 + 1))
    normalization = (2.0 * l_value + 1.0) / (2.0 * j_value + 1.0)
    norm2 = 0.0
    for m1 in helicities1:
        for m2 in helicities2:
            helicity = m1 - m2
            if abs(helicity) > j_value:
                continue
            coefficient = clebsch_gordan(
                float(l_value), s_value, 0.0, float(helicity), j_value, float(helicity)
            ) * clebsch_gordan(
                float(j1), float(j2), float(m1), -float(m2), s_value, float(helicity)
            )
            norm2 += normalization * coefficient * coefficient
    tol = 1.0e-12
    return norm2 > tol * tol


# Compute normalized spin-rank polarization in the occupation basis
@functools.cache
def _occupations(rank: int, projection: int) -> dict[tuple[int, int, int], float]:
    state: dict[tuple[int, int, int], float] = {(rank, 0, 0): 1.0}
    for current in range(rank, projection, -1):
        denominator = math.sqrt(float((rank + current) * (rank - current + 1)))
        lowered: dict[tuple[int, int, int], float] = {}
        for (plus, zero, minus), coefficient in state.items():
            if plus > 0:
                key = (plus - 1, zero + 1, minus)
                lowered[key] = (
                    lowered.get(key, 0.0)
                    + coefficient * math.sqrt(2.0 * float(plus * (zero + 1))) / denominator
                )
            if zero > 0:
                key = (plus, zero - 1, minus + 1)
                lowered[key] = (
                    lowered.get(key, 0.0)
                    + coefficient * math.sqrt(2.0 * float(zero * (minus + 1))) / denominator
                )
        state = lowered
    return state


# Compute ordered Cartesian component of a normalized symmetric tensor
def _tensor_component(
    rank: int, projection: int, occupation: tuple[int, int, int]
) -> float:
    coefficient = _occupations(rank, projection).get(occupation, 0.0)
    if coefficient == 0.0:
        return 0.0
    log_weight = 0.5 * (
        sum(math.lgamma(float(value + 1)) for value in occupation) - math.lgamma(float(rank + 1))
    )
    return coefficient * math.exp(log_weight)


# Compute the raw Cartesian STF coupling relative to normalized SU(2) CG
def stf_norm(rank1: int, rank2: int, output_rank: int) -> float:
    if output_rank > rank1 + rank2 or output_rank < abs(rank1 - rank2):
        raise ValueError("STF ranks violate the triangle rule")
    difference = rank1 + rank2 - output_rank
    contractions = difference // 2
    epsilon = difference % 2
    projection = output_rank - rank1
    occupation = (rank2 - contractions - epsilon, epsilon, contractions)
    cartesian = _tensor_component(rank2, projection, occupation)
    cg = clebsch_gordan(
        float(rank1),
        float(rank2),
        float(rank1),
        float(projection),
        float(output_rank),
        float(output_rank),
    )
    if abs(cartesian) <= 1.0e-12 or abs(cg) <= 1.0e-12:
        raise ValueError("STF highest-weight projection is zero")
    return abs(cartesian / cg)


# Compute the raw pole normalization multiplying one normalized JW LS tensor
def ls_norm(
    mother: Particle, p1: Particle, p2: Particle, l_value: int, s_value: float
) -> float:
    if mother.spin2 < 0 or mother.spin2 % 2 != 0 or not is_integer(s_value):
        raise ValueError("raw pole normalization requires integer mother spin and S")
    j_value = mother.spin2 // 2
    s_int = int(round(s_value))
    if p1.spin2 % 2 == 0 and p2.spin2 % 2 == 0:
        leg_normalization = stf_norm(p1.spin2 // 2, p2.spin2 // 2, s_int)
    elif p1.spin2 == 1 and p2.spin2 == 1:
        leg_normalization = math.sqrt(2.0)
    else:
        raise ValueError("raw pole normalization does not support mixed spinor and tensor legs")
    orbital = math.exp(
        0.5
        * (
            float(l_value) * math.log(2.0)
            + 2.0 * math.lgamma(float(l_value + 1))
            - math.lgamma(float(2 * l_value + 1))
        )
    )
    return (
        leg_normalization
        * stf_norm(l_value, s_int, j_value)
        * orbital
        * math.sqrt((2.0 * j_value + 1.0) / (2.0 * l_value + 1.0))
    )


# Compute compact independent physical-pole helicity rows for one LS reference
def helicity_rows(
    mother: Particle,
    p1: Particle,
    p2: Particle,
    allowed: list[tuple[int, float]],
    reference: tuple[int, float],
    tolerance: float = 1.0e-12,
    *,
    p_symmetry: bool = True,
) -> list[list[Any]]:
    if reference not in allowed:
        raise ValueError("canonical helicity reference is not allowed")
    j_value = 0.5 * mother.spin2
    j1 = 0.5 * p1.spin2
    j2 = 0.5 * p2.spin2
    helicities1 = (
        (-1.0, 1.0) if p1.pdg == 22 else tuple(-j1 + float(index) for index in range(p1.spin2 + 1))
    )
    helicities2 = (
        (-1.0, 1.0) if p2.pdg == 22 else tuple(-j2 + float(index) for index in range(p2.spin2 + 1))
    )
    raw = ls_norm(mother, p1, p2, reference[0], reference[1])
    coordinates: list[tuple[float, float, complex]] = []
    for m1 in helicities1:
        for m2 in helicities2:
            helicity = m1 - m2
            support = False
            for l_value, s_value in allowed:
                coefficient = clebsch_gordan(
                    float(l_value), s_value, 0.0, float(helicity), j_value, float(helicity)
                ) * clebsch_gordan(j1, j2, float(m1), -float(m2), s_value, float(helicity))
                support = support or abs(coefficient) > tolerance
            if not support:
                continue
            coefficient = raw * math.sqrt((2.0 * reference[0] + 1.0) / (2.0 * j_value + 1.0))
            coefficient *= clebsch_gordan(
                float(reference[0]),
                reference[1],
                0.0,
                float(helicity),
                j_value,
                float(helicity),
            )
            coefficient *= clebsch_gordan(
                j1,
                j2,
                float(m1),
                -float(m2),
                reference[1],
                float(helicity),
            )
            coordinates.append((m1, m2, complex(coefficient, 0.0)))

    rows: list[list[Any]] = []
    covered: set[tuple[float, float]] = set()
    for m1, m2, coefficient in coordinates:
        if (m1, m2) in covered:
            continue
        magnitude, phase = mag_phase(coefficient.real, coefficient.imag, tolerance)
        rows.append(
            [
                clean_number(m1),
                clean_number(m2),
                clean_number(magnitude),
                clean_number(phase),
            ]
        )
        covered.add((m1, m2))
        if p_symmetry:
            covered.add((-m1, -m2))
    return rows


# Compute true if one fixed spin LS row is allowed
def ls_allowed(
    mother: Particle,
    p1: Particle,
    p2: Particle,
    l_value: int,
    s_value: float,
    c_symmetry: bool,
    p_symmetry: bool,
    *,
    crossed: bool = False,
) -> bool:
    j_value = 0.5 * mother.spin2
    s1 = 0.5 * p1.spin2
    s2 = 0.5 * p2.spin2
    if not is_integer(s1 + s2 + s_value):
        return False
    if s_value < abs(s1 - s2) - 1e-9 or s_value > s1 + s2 + 1e-9:
        return False
    if j_value < abs(float(l_value) - s_value) - 1e-9 or j_value > float(l_value) + s_value + 1e-9:
        return False
    if p_symmetry and p1.parity * p2.parity * (1 if l_value % 2 == 0 else -1) != mother.parity:
        return False
    if not c_allowed(mother, p1, p2, l_value, s_value, c_symmetry and not crossed):
        return False
    if p1.pdg == p2.pdg and (not crossed or (has_cparity(p1) and has_cparity(p2))):
        if not is_integer(s_value):
            return False
        s_int = int(round(s_value))
        if p1.spin2 % 2 == 0:
            if (l_value - s_int) % 2 != 0:
                return False
        elif (l_value + s_int) % 2 != 0:
            return False
    return ls_supported(mother, p1, p2, l_value, s_value)


# Compute allowed physical pole LS rows
def ls_rows(
    mother: Particle,
    p1: Particle,
    p2: Particle,
    c_symmetry: bool,
    p_symmetry: bool,
    *,
    crossed: bool = False,
) -> list[tuple[int, float]]:
    out: list[tuple[int, float]] = []
    max_s2 = max(0, p1.spin2 + p2.spin2)
    for two_s in range(max_s2 + 1):
        s_value = 0.5 * two_s
        max_l = int(math.floor(0.5 * mother.spin2 + s_value + 1e-9))
        for l_value in range(max_l + 1):
            allowed = ls_allowed(
                mother,
                p1,
                p2,
                l_value,
                s_value,
                c_symmetry,
                p_symmetry,
                crossed=crossed,
            )
            if allowed:
                out.append((l_value, s_value))
    return out
