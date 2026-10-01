# Tool functions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
import re
from fractions import Fraction
from functools import cache

from core.numerics import array

PHASE_U_SUFFIX = "@PHASE_U"
PHASE_V_SUFFIX = "@PHASE_V"

COMPLEX_RE_SUFFIX = "@RE"
COMPLEX_IM_SUFFIX = "@IM"
MAGNITUDE_SUFFIX = "@MAG"
COMPLEX_EPS = 1.0e-12

PROJECTIVE_VECTOR_SUFFIX = "@PROJECTIVE"
SPHERICAL_VECTOR_SUFFIX = "@SPHERICAL"
PROJECTIVE_NORM_SUFFIX = "@NORM"
PROJECTIVE_COEFF_RE_SUFFIX = "@COEFF_RE"
PROJECTIVE_COEFF_IM_SUFFIX = "@COEFF_IM"
AJZ_PROJECTIVE_SUFFIX = "@AJZP"
SIMPLEX_VECTOR_SUFFIX = "@SIMPLEX"
RAW_PHASE_SUFFIX = "@PHASE"
PROJECTIVE_ANGLE_LIMIT = 0.5 * math.pi

TUNE_MODES = {
    "none",
    "magnitude",
    "phase_raw",
    "phase_cayley",
    "mag_phase_raw",
    "mag_phase_cartesian",
    "mag_phase_cayley",
}


# String handling


def find_str_between(s: str, start: str, end: str) -> str:
    i = s.find(start)
    if i < 0:
        raise ValueError(f"find_str_between: start delimiter '{start}' not found")
    i += len(start)

    j = s.find(end, i)
    if j < 0:
        raise ValueError(f"find_str_between: end delimiter '{end}' not found")

    return s[i:j]


@cache
def _compiled_split_pattern(delimiters: tuple[str, ...]):
    """Return a compiled regex for the given delimiters."""
    return re.compile("|".join(map(re.escape, delimiters)))


def substring_split(s: str, delimiters: list[str], maxsplit: int = 0) -> list[str]:
    """
    Split with multiple delimiters

    Example:
        substring_split('DECAY|1234:[321,-321]:zeta.GP', ['|',':'])
        # Output: ['DECAY', '1234', '[321,-321]', 'zeta.GP']
    """
    if len(delimiters) == 0:
        raise ValueError("substring_split: delimiters must be non-empty")
    return _compiled_split_pattern(tuple(delimiters)).split(s, maxsplit)


def string_to_index(s: str, delim: str) -> tuple[int]:
    """
    Returns tuple with indices given a delimiter
    """
    parts = s.split(delim)
    if any(item == "" for item in parts):
        raise ValueError(f"string_to_index: malformed index string '{s}'")
    try:
        return tuple(int(item) for item in parts)
    except ValueError as e:
        raise ValueError(f"string_to_index: non-integer token in '{s}'") from e


# Phase and magnitude/phase optimizer encodings


# Wrap a real phase with the same branch convention for NumPy and Torch
def canonicalize_phase(phi: float) -> float:
    phi = array.asarray(phi, dtype=float)
    xp = array.namespace(phi)
    value = xp.arctan2(xp.sin(phi), xp.cos(phi))
    return value - 2.0 * math.pi if value >= math.pi else value


# Compute the canonical phase of one signed real amplitude
def real_phase(value: float) -> float:
    coefficient = array.asarray(value, dtype=float)
    xp = array.namespace(coefficient)
    if not xp.isfinite(coefficient):
        raise ValueError("real_phase: coefficient must be finite")
    return xp.where(coefficient >= 0.0, array.asarray(0.0, like=coefficient), array.asarray(-math.pi, like=coefficient))[()]


def encode_phase(phi: float) -> tuple[float, float]:
    """
    Return symmetric Cayley variables u,v for a scalar phase
    """
    value = canonicalize_phase(phi)
    t = math.tan(0.25 * value)
    return t, t


def decode_phase_components(u: float, v: float) -> float:
    """
    Decode two Cayley phase variables into a canonical scalar phase

    z(u,v) = ((1+iu)/(1-iu))*((1+iv)/(1-iv))
           = (1-uv + i(u+v)) / (1-uv-i(u+v))
    """
    u = array.asarray(u, like=v, dtype=float)
    v = array.asarray(v, like=u)
    a = 1.0 - u * v
    b = u + v

    return canonicalize_phase(2.0 * array.namespace(u).arctan2(b, a))


def phase_u_key(base_key: str) -> str:
    return f"{base_key}{PHASE_U_SUFFIX}"


def phase_v_key(base_key: str) -> str:
    return f"{base_key}{PHASE_V_SUFFIX}"


# Compute whether a key has one of the selected suffixes
def _has_suffix(key: str, suffixes: tuple[str, ...]) -> bool:
    return str(key).endswith(suffixes)


# Remove the matching encoded suffix from a key
def _base_before_suffix(key: str, suffixes: tuple[str, ...], label: str) -> str:
    key = str(key)
    suffix = next((item for item in suffixes if key.endswith(item)), None)
    if suffix is None:
        raise ValueError(f'{label}: key "{key}" has no recognized suffix')
    return key[: -len(suffix)]


def is_phase_encoded_key(key: str) -> bool:
    return _has_suffix(key, (PHASE_U_SUFFIX, PHASE_V_SUFFIX))


def encoded_phase_base_key(key: str) -> str:
    return _base_before_suffix(key, (PHASE_U_SUFFIX, PHASE_V_SUFFIX), "encoded_phase_base_key")


def encode_polar(mag: float, phi: float) -> tuple[float, float]:
    """Encode a polar complex number as Cartesian components."""
    value = canonicalize_phase(float(phi))
    radius = float(mag)
    return radius * math.cos(value), radius * math.sin(value)


# Decode Cartesian amplitudes with finite derivatives at a vanishing coefficient
def decode_cartesian(u: float, v: float, *, eps: float = COMPLEX_EPS) -> tuple[float, float]:
    u = array.asarray(u, like=v, dtype=float)
    v = array.asarray(v, like=u)
    xp = array.namespace(u)
    magnitude = array.hypot(u, v)
    active = magnitude >= eps
    phase = canonicalize_phase(xp.arctan2(xp.where(active, v, 0.0), xp.where(active, u, 1.0)))
    return magnitude if active else magnitude * 0.0, phase


def complex_re_key(base_key: str) -> str:
    return f"{base_key}{COMPLEX_RE_SUFFIX}"


def complex_im_key(base_key: str) -> str:
    return f"{base_key}{COMPLEX_IM_SUFFIX}"


# Compute one explicitly marked physical magnitude key
def magnitude_key(base_key: str) -> str:
    return f"{base_key}{MAGNITUDE_SUFFIX}"


# Compute true for explicitly marked physical magnitude keys
def is_magnitude_key(key: str) -> bool:
    return str(key).endswith(MAGNITUDE_SUFFIX)


# Compute the coupling name of one explicitly marked magnitude key
def magnitude_base_key(key: str) -> str:
    return _base_before_suffix(key, (MAGNITUDE_SUFFIX,), "magnitude_base_key")


def is_complex_encoded_key(key: str) -> bool:
    return _has_suffix(key, (COMPLEX_RE_SUFFIX, COMPLEX_IM_SUFFIX))


def encoded_complex_base_key(key: str) -> str:
    return _base_before_suffix(key, (COMPLEX_RE_SUFFIX, COMPLEX_IM_SUFFIX), "encoded_complex_base_key")


def validate_tune_mode(mode: str, *, name: str) -> str:
    mode = str(mode).strip()
    if mode not in TUNE_MODES:
        raise ValueError(f'Unknown {name} "{mode}"; expected one of {sorted(TUNE_MODES)}')
    return mode


# Compute one indexed optimizer parameter key
def _indexed_key(base_key: str, angle_index: int, marker: str, suffix: str) -> str:
    return f"{base_key}{marker}{int(angle_index)}{suffix}"


# Parse one indexed optimizer parameter key
def _parse_indexed_key(
    key: str, *, marker: str, suffix: str, label: str, required: str | None = None
) -> tuple[str, int]:
    key = str(key)
    if not key.endswith(suffix) or (required is not None and required not in key):
        raise ValueError(f'{label}: key "{key}" is not encoded with {suffix}')
    base_key, found, index = key[: -len(suffix)].rpartition(marker)
    if not found or not base_key:
        raise ValueError(f'{label}: malformed indexed key "{key}"')
    try:
        angle_index = int(index)
    except ValueError as exc:
        raise ValueError(f'{label}: non-integer index in "{key}"') from exc
    if angle_index < 0:
        raise ValueError(f'{label}: negative index in "{key}"')
    return base_key, angle_index


# Compute a real-projective hyperspherical angle parameter key
def projective_angle_key(base_key: str, angle_index: int) -> str:
    return _indexed_key(base_key, angle_index, "_angle", PROJECTIVE_VECTOR_SUFFIX)


# Compute true for real-projective hyperspherical angle parameter keys
def is_projective_angle_key(key: str) -> bool:
    return str(key).endswith(PROJECTIVE_VECTOR_SUFFIX)


# Parse a real-projective hyperspherical angle parameter key
def parse_projective_angle_key(key: str) -> tuple[str, int | str]:
    key = str(key)
    stem = key[: -len(PROJECTIVE_VECTOR_SUFFIX)] if is_projective_angle_key(key) else ""
    physical = re.fullmatch(r"(.*:(?:g_ls|alpha_ls|helicity|g_tensor))\((.*)\)", stem)
    if physical is not None:
        return physical.group(1), physical.group(2)
    return _parse_indexed_key(key, marker="_angle", suffix=PROJECTIVE_VECTOR_SUFFIX, label="parse_projective_angle_key")


# Compute a real spherical hyperspherical angle parameter key
def spherical_angle_key(base_key: str, angle_index: int) -> str:
    return _indexed_key(base_key, angle_index, "_angle", SPHERICAL_VECTOR_SUFFIX)


# Compute true for real spherical hyperspherical angle parameter keys
def is_spherical_angle_key(key: str) -> bool:
    return str(key).endswith(SPHERICAL_VECTOR_SUFFIX)


# Parse a real spherical hyperspherical angle parameter key
def parse_spherical_angle_key(key: str) -> tuple[str, int | str]:
    key = str(key)
    stem = key[: -len(SPHERICAL_VECTOR_SUFFIX)] if is_spherical_angle_key(key) else ""
    physical = re.fullmatch(r"(.*:(?:g_ls|alpha_ls|helicity|g_tensor))\((.*)\)", stem)
    if physical is not None:
        return physical.group(1), physical.group(2)
    return _parse_indexed_key(key, marker="_angle", suffix=SPHERICAL_VECTOR_SUFFIX, label="parse_spherical_angle_key")


# Compute an overall norm key carrying the projective reference row label
def projective_norm_key(reference_key: str) -> str:
    return f"{reference_key}{PROJECTIVE_NORM_SUFFIX}"


# Compute true for an overall projective-vector norm key
def is_projective_norm_key(key: str) -> bool:
    return str(key).endswith(PROJECTIVE_NORM_SUFFIX)


# Parse an overall projective norm key into its vector base and reference label
def parse_projective_norm_key(key: str) -> tuple[str, str]:
    key = str(key)
    stem = key[: -len(PROJECTIVE_NORM_SUFFIX)] if is_projective_norm_key(key) else ""
    physical = re.fullmatch(r"(.*:(?:g_ls|alpha_ls|helicity|g_tensor))\((.*)\)", stem)
    if physical is None:
        raise ValueError(f'parse_projective_norm_key: malformed norm key "{key}"')
    return physical.group(1), physical.group(2)


# Compute Cartesian keys for one complex coefficient multiplying a projective vector
def projective_coefficient_keys(reference_key: str) -> tuple[str, str]:
    return f"{reference_key}{PROJECTIVE_COEFF_RE_SUFFIX}", f"{reference_key}{PROJECTIVE_COEFF_IM_SUFFIX}"


# Compute true for Cartesian projective-vector coefficient keys
def is_projective_coefficient_key(key: str) -> bool:
    return str(key).endswith((PROJECTIVE_COEFF_RE_SUFFIX, PROJECTIVE_COEFF_IM_SUFFIX))


# Parse a Cartesian projective-vector coefficient key
def parse_projective_coefficient_key(key: str) -> tuple[str, str, str]:
    key = str(key)
    suffix = next(
        (item for item in (PROJECTIVE_COEFF_RE_SUFFIX, PROJECTIVE_COEFF_IM_SUFFIX) if key.endswith(item)),
        None,
    )
    stem = key[: -len(suffix)] if suffix is not None else ""
    physical = re.fullmatch(r"(.*:(?:g_ls|alpha_ls|helicity|g_tensor))\((.*)\)", stem)
    if physical is None:
        raise ValueError(f'parse_projective_coefficient_key: malformed Cartesian coefficient key "{key}"')
    component = "re" if suffix == PROJECTIVE_COEFF_RE_SUFFIX else "im"
    return physical.group(1), physical.group(2), component


# Compute the common physical phase associated with one signed-real direction
def projective_phase_base_key(base_key: str) -> str | None:
    parts = substring_split(s=str(base_key), delimiters=["|", ":"])
    if (
        len(parts) == 4
        and parts[0] == "RES"
        and parts[2] in {"MP", "XP", "GP", "TP"}
        and parts[3] in {"g_ls", "helicity", "g_tensor"}
    ):
        return f"RES|{parts[1]}:{parts[2]}:phi"
    # Shared decay shapes have independent phases for each Pomeron model
    return None


# Compute continuum symmetry multiplicities from physical compact row labels
def projective_completion_weights(base_key: str, reference: str, angle_names: list[str]) -> list[float] | None:
    labels = [str(reference), *[str(parse_projective_angle_key(name)[1]) for name in angle_names]]
    return continuum_completion_weights(base_key, labels)


# Parse one direct GP resonance row into its vector base and physical row label
def gp_res_component_key(key: str, suffixes: tuple[str, ...]) -> tuple[str, str] | None:
    stem = _base_before_suffix(key, suffixes, "gp_res_component_key")
    physical = re.fullmatch(r"(RES\|[^:]+:GP:(?:g_ls|helicity))\((.*)\)", stem)
    return None if physical is None else (physical.group(1), physical.group(2))


# Count the C++ continuum parity and identical-leg symmetry orbits
def continuum_completion_weights(base_key: str, labels: list[str]) -> list[float] | None:
    parts = substring_split(s=str(base_key), delimiters=["|", ":"])
    if len(parts) != 4 or parts[0] not in {"CON_MP", "CON_XP", "CON_GP"} or parts[3] not in {"helicity", "g_ls"}:
        return None
    pair, _, sector = parts[2].partition("/")
    pdgs = [int(pdg) for pdg in pair.strip("[]").split(",")]
    identical = sector == "self" and pdgs[0] == pdgs[1]
    weights = []
    for label in labels:
        row = tuple(Fraction(value.strip()) for value in label.split(","))
        if parts[3] == "g_ls":
            weights.append(2.0 if len(row) == 3 and row[2] else 1.0)
        else:
            orbit = {row, tuple(-value for value in row)}
            if identical:
                orbit.update((item[1], item[0], *item[2:]) for item in tuple(orbit))
            weights.append(float(len(orbit)))
    return weights


# Normalize a real vector after choosing one representative of its antipodal pair
def canonicalize_projective_vector(values: list[float], *, eps: float = COMPLEX_EPS) -> list[float]:
    vector = normalize_spherical_vector(values, eps=eps)
    sign_reference = vector[0]
    if abs(sign_reference) <= eps:
        sign_reference = next((value for value in vector[1:] if abs(value) > eps), 0.0)
    if sign_reference < 0.0:
        vector = [-value for value in vector]
    return vector


# Decode hyperspherical coordinates into an oriented unit vector
def _hyperspherical_vector(angles, *, sphere=False):
    angles = array.stack(list(angles))
    xp = array.namespace(angles)
    values = [0.0] * (len(angles) + 1)
    remaining = 1.0
    for k, angle in enumerate(angles):
        limit = math.pi if sphere and k + 1 == len(angles) else PROJECTIVE_ANGLE_LIMIT
        if not xp.isfinite(angle) or abs(angle) > limit + COMPLEX_EPS:
            raise ValueError(f"Hyperspherical angle {k} outside [-{limit}, {limit}]")
        values[k + 1] = remaining * xp.sin(angle)
        remaining = remaining * xp.cos(angle)
    values[0] = remaining
    return normalize_spherical_vector(values)


# Decode signed hyperspherical angles into a normalized projective vector
def projective_vector_from_angles(angles: list[float]) -> list[float]:
    return canonicalize_projective_vector(_hyperspherical_vector(angles))


# Decode signed hyperspherical angles without changing the hemisphere orientation
def hemisphere_vector_from_angles(angles: list[float]) -> list[float]:
    return _hyperspherical_vector(angles)


# Encode a real projective vector as signed hyperspherical angles
def projective_angles_from_vector(values: list[float]) -> list[float]:
    v = canonicalize_projective_vector(values)
    angles = []
    for k in range(len(v) - 1):
        remaining = math.sqrt(v[0] * v[0] + sum(value * value for value in v[k + 1 :]))
        if remaining <= COMPLEX_EPS:
            angles.append(0.0)
            continue
        argument = max(-1.0, min(1.0, v[k + 1] / remaining))
        angles.append(math.asin(argument))
    return angles


# Normalize a real vector while retaining its global orientation
def normalize_spherical_vector(values: list[float], *, eps: float = COMPLEX_EPS) -> list[float]:
    vector = array.stack(list(values))
    xp = array.namespace(vector)
    if not len(vector) or not xp.all(xp.isfinite(vector)):
        raise ValueError("normalize_spherical_vector: expected a finite non-empty vector")
    norm = xp.linalg.norm(vector)
    if norm <= eps:
        raise ValueError("normalize_spherical_vector: zero vector")
    return list(vector / norm)


# Decode full sphere hyperspherical angles into a normalized real vector
def spherical_vector_from_angles(angles: list[float]) -> list[float]:
    return _hyperspherical_vector(angles, sphere=True)


# Encode an oriented real unit vector as full sphere hyperspherical angles
def spherical_angles_from_vector(values: list[float]) -> list[float]:
    vector = normalize_spherical_vector(values)
    if len(vector) == 1:
        return []
    angles = []
    for k in range(len(vector) - 2):
        remaining = math.sqrt(vector[0] * vector[0] + sum(value * value for value in vector[k + 1 :]))
        if remaining <= COMPLEX_EPS:
            angles.append(0.0)
            continue
        argument = max(-1.0, min(1.0, vector[k + 1] / remaining))
        angles.append(math.asin(argument))
    angles.append(canonicalize_phase(math.atan2(vector[-1], vector[0])))
    return angles


# Decode spherical angles directly into oriented unit vector features
def spherical_embedding_from_angles(angles: list[float]) -> list[float]:
    return spherical_vector_from_angles(angles)


# Compute stable labels for oriented spherical vector components
def spherical_embedding_feature_names(base_key: str, vector_size: int) -> list[str]:
    if int(vector_size) < 1:
        raise ValueError("spherical_embedding_feature_names: vector size must be positive")
    return [f"{base_key}[{index}]" for index in range(int(vector_size))]


# Embed a normalized direction through the sign-invariant projector vv^T
def projective_embedding_from_vector(values: list[float]) -> list[float]:
    vector = canonicalize_projective_vector(values)
    features = []
    for i, left in enumerate(vector):
        for j in range(i, len(vector)):
            scale = 1.0 if i == j else math.sqrt(2.0)
            features.append(scale * left * vector[j])
    return features


# Decode projective angles directly into sign-invariant projector features
def projective_embedding_from_angles(angles: list[float]) -> list[float]:
    return projective_embedding_from_vector(projective_vector_from_angles(angles))


# Compute stable labels for the independent symmetric-projector components
def projective_embedding_feature_names(base_key: str, vector_size: int) -> list[str]:
    vector_size = int(vector_size)
    if vector_size < 1:
        raise ValueError("projective_embedding_feature_names: vector size must be positive")
    return [f"{base_key}[{i},{j}]" for i in range(vector_size) for j in range(i, vector_size)]


# Compute a coherent MP hyperspherical angle key
def ajzp_angle_key(base_key: str, angle_index: int) -> str:
    return _indexed_key(base_key, angle_index, "_angle", AJZ_PROJECTIVE_SUFFIX)


# Compute true for coherent MP angle keys
def is_ajzp_angle_key(key: str) -> bool:
    return str(key).endswith(AJZ_PROJECTIVE_SUFFIX) and ".a_Jz_angle" in str(key)


# Parse a coherent MP hyperspherical angle key
def parse_ajzp_angle_key(key: str) -> tuple[str, int]:
    return _parse_indexed_key(
        key, marker="_angle", suffix=AJZ_PROJECTIVE_SUFFIX, label="parse_ajzp_angle_key", required=".a_Jz_angle"
    )


# Compute a probability-simplex angle parameter key
def simplex_theta_key(base_key: str, angle_index: int) -> str:
    return _indexed_key(base_key, angle_index, "_theta", SIMPLEX_VECTOR_SUFFIX)


# Compute true for probability-simplex angle parameter keys
def is_simplex_theta_key(key: str) -> bool:
    return str(key).endswith(SIMPLEX_VECTOR_SUFFIX)


# Parse a probability-simplex angle parameter key
def parse_simplex_theta_key(key: str) -> tuple[str, int]:
    return _parse_indexed_key(key, marker="_theta", suffix=SIMPLEX_VECTOR_SUFFIX, label="parse_simplex_theta_key")


# Decode first-orthant angles into normalized probability populations
def simplex_probabilities_from_angles(angles: list[float]) -> list[float]:
    amplitudes = []
    sinprod = 1.0
    for angle in angles:
        angle = array.asarray(angle, dtype=float)
        xp = array.namespace(angle)
        if angle < 0.0 or angle > 0.5 * math.pi:
            raise ValueError("simplex_probabilities_from_angles: angle outside [0, pi/2]")
        amplitudes.append(sinprod * xp.cos(angle))
        sinprod = sinprod * xp.sin(angle)
    amplitudes.append(sinprod)
    return [value * value for value in amplitudes]


# Encode normalized probability populations into first-orthant angles
def simplex_angles_from_probabilities(probabilities: list[float]) -> list[float]:
    values = [float(value) for value in probabilities]
    if any(value < 0.0 for value in values):
        raise ValueError("simplex_angles_from_probabilities: negative probability")
    total = sum(values)
    if total <= 0.0:
        raise ValueError("simplex_angles_from_probabilities: zero total probability")
    amplitudes = [math.sqrt(value / total) for value in values]
    angles = []
    for k in range(len(amplitudes) - 1):
        radius = math.sqrt(sum(value * value for value in amplitudes[k:]))
        arg = 1.0 if radius <= 0.0 else amplitudes[k] / radius
        angles.append(math.acos(max(-1.0, min(1.0, arg))))
    return angles


# Compute a marked optimizer key for one physical phase
def raw_phase_key(base_key: str) -> str:
    return f"{base_key}{RAW_PHASE_SUFFIX}"


# Compute true for marked physical phase optimizer keys
def is_raw_phase_key(key: str) -> bool:
    return str(key).endswith(RAW_PHASE_SUFFIX)


# Compute the steering-card target of one marked phase key
def raw_phase_base_key(key: str) -> str:
    key = str(key)
    if not is_raw_phase_key(key):
        raise ValueError(f'raw_phase_base_key: key "{key}" is not a marked phase')
    return key[: -len(RAW_PHASE_SUFFIX)]


# Build optional topology groups from marked optimizer parameter names
def build_parameter_topology(parameter_names) -> dict:
    names = [str(name) for name in parameter_names]
    groups = []
    grouped = {}
    phase_names = {}
    norm_names = {}
    cartesian_coefficients = {}
    magnitude_components = {}
    for name in names:
        if is_projective_angle_key(name):
            base, index = parse_projective_angle_key(name)
            grouped.setdefault(("projective", base), []).append((index, name))
        elif is_spherical_angle_key(name):
            base, index = parse_spherical_angle_key(name)
            grouped.setdefault(("sphere", base), []).append((index, name))
        elif is_ajzp_angle_key(name):
            base, index = parse_ajzp_angle_key(name)
            grouped.setdefault(("projective", base), []).append((index, name))
        elif is_simplex_theta_key(name):
            base, index = parse_simplex_theta_key(name)
            grouped.setdefault(("simplex", base), []).append((index, name))
        elif is_raw_phase_key(name):
            phase_names[raw_phase_base_key(name)] = name
        elif is_projective_norm_key(name):
            base, _ = parse_projective_norm_key(name)
            if base in norm_names:
                raise ValueError(f'build_parameter_topology: duplicate projective norm for "{base}"')
            norm_names[base] = name
        elif is_projective_coefficient_key(name):
            base, _, component = parse_projective_coefficient_key(name)
            if component in cartesian_coefficients.setdefault(base, {}):
                raise ValueError(f'build_parameter_topology: duplicate Cartesian projective coefficient for "{base}"')
            cartesian_coefficients[base][component] = name
        elif is_magnitude_key(name):
            parsed = gp_res_component_key(name, (MAGNITUDE_SUFFIX,))
            if parsed is not None:
                base, label = parsed
                magnitude_components.setdefault(base, []).append((label, name))

    consumed_phases = set()
    direction_groups = {}
    for (kind, base), entries in sorted(grouped.items()):
        numeric = all(isinstance(index, int) for index, _ in entries)
        physical = all(isinstance(index, str) for index, _ in entries)
        if not numeric and not physical:
            raise ValueError(f'build_parameter_topology: mixed labels in {kind} group "{base}"')
        if numeric:
            sorted_entries = sorted(entries)
            expected = list(range(len(sorted_entries)))
            found = [index for index, _ in sorted_entries]
            if found != expected:
                raise ValueError(f'build_parameter_topology: incomplete {kind} group "{base}"')
        else:
            sorted_entries = sorted(entries, key=lambda entry: names.index(entry[1]))
            found = [index for index, _ in sorted_entries]
            if len(found) != len(set(found)):
                raise ValueError(f'build_parameter_topology: duplicate labels in {kind} group "{base}"')
        ordered = [name for _, name in sorted_entries]
        direction_groups[(kind, base)] = ordered

    consumed_norms = set()
    spherical_bases = {base for kind, base in direction_groups if kind == "sphere"}
    for base in sorted(spherical_bases):
        angles = direction_groups.get(("sphere", base), [])
        phase_base = projective_phase_base_key(base)
        phase_name = phase_names.get(phase_base)
        norm_name = norm_names.get(base)
        if norm_name is None or phase_name is None:
            continue
        direction_groups.pop(("sphere", base))
        consumed_phases.add(phase_base)
        consumed_norms.add(base)
        groups.append(
            {
                "kind": "polar_sphere",
                "base": base,
                "parameters": [norm_name, phase_name, *angles],
                "period": 2.0 * math.pi,
            }
        )

    projective_bases = (
        {base for kind, base in direction_groups if kind == "projective"}
        | (set(norm_names) - consumed_norms)
        | set(cartesian_coefficients)
    )
    for base in sorted(projective_bases):
        angles = direction_groups.pop(("projective", base), [])
        phase_base = projective_phase_base_key(base)
        phase_name = phase_names.get(phase_base)
        cartesian = cartesian_coefficients.get(base)
        norm_name = norm_names.get(base)
        if cartesian is not None:
            if set(cartesian) != {"re", "im"}:
                raise ValueError(f'build_parameter_topology: incomplete Cartesian projective coefficient for "{base}"')
            if norm_name is not None or phase_name is not None:
                raise ValueError(f'build_parameter_topology: mixed projective residue coordinates for "{base}"')
            _, reference, _ = parse_projective_coefficient_key(cartesian["re"])
            group = {
                "kind": "cartesian_projective",
                "base": base,
                "parameters": [cartesian["re"], cartesian["im"], *angles],
            }
            if (weights := projective_completion_weights(base, reference, angles)) is not None:
                group["weights"] = weights
            groups.append(group)
            continue
        if norm_name is not None and base.startswith(("CON_MP|", "CON_XP|", "CON_GP|")):
            # A single helicity coupling has only a scalar norm
            if not angles:
                continue
            _, reference = parse_projective_norm_key(norm_name)
            group = {
                "kind": "radial_projective",
                "base": base,
                "parameters": [norm_name, *angles],
            }
            if (weights := projective_completion_weights(base, reference, angles)) is not None:
                group["weights"] = weights
            groups.append(group)
            continue
        if norm_name is not None and phase_name is not None:
            consumed_phases.add(phase_base)
            groups.append(
                {
                    "kind": "polar_projective",
                    "base": base,
                    "parameters": [norm_name, phase_name, *angles],
                    "period": 2.0 * math.pi,
                }
            )
            continue
        if phase_name is not None and angles:
            consumed_phases.add(phase_base)
            groups.append(
                {
                    "kind": "phase_projective",
                    "base": base,
                    "parameters": [phase_name, *angles],
                    "period": 2.0 * math.pi,
                }
            )
            continue
        if angles:
            group = {"kind": "projective", "base": base, "parameters": angles}
            if norm_name is not None:
                _, reference = parse_projective_norm_key(norm_name)
                if (weights := projective_completion_weights(base, reference, angles)) is not None:
                    group["weights"] = weights
            groups.append(group)

    for base, entries in magnitude_components.items():
        ordered = sorted(entries, key=lambda entry: names.index(entry[1]))
        parameters = [entry[1] for entry in ordered]
        if base.startswith("RES|"):
            phase_base = projective_phase_base_key(base)
            phase_name = phase_names.get(phase_base)
            if phase_name is None:
                continue
            consumed_phases.add(phase_base)
            groups.append(
                {
                    "kind": "polar_components",
                    "base": base,
                    "parameters": [phase_name, *parameters],
                    "period": 2.0 * math.pi,
                }
            )

    for (kind, base), ordered in direction_groups.items():
        groups.append({"kind": kind, "base": base, "parameters": ordered})
    for base, name in phase_names.items():
        if base not in consumed_phases:
            groups.append({"kind": "phase", "base": base, "parameters": [name], "period": 2.0 * math.pi})
    groups.sort(key=lambda group: min(names.index(name) for name in group["parameters"]))
    return {"schema_version": 1, "groups": groups} if groups else {}


def tune_mode_has_magnitude(mode: str) -> bool:
    return mode in {"magnitude", "mag_phase_raw", "mag_phase_cayley"}


def tune_mode_is_cartesian(mode: str) -> bool:
    return mode == "mag_phase_cartesian"
