# Read and project coupled channel eikonal matrix outputs
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from core.io.serialize import decode_doubles, load_json_file

_CACHE_VERSION = 1
_SPIN_ENTRIES = 16
_SPIN_DIAGONAL = (0, 5, 10, 15)
_FILENAME = re.compile(r"^MATRIX_N(?P<n>\d+)_(?P<b1>-?\d+)_(?P<b2>-?\d+)_(?P<hash>\d+)\.json$")


@dataclass(frozen=True)
class EikonalMatrixFile:
    path: Path
    nchannels: int
    beam1: int
    beam2: int
    sqrts: float


@dataclass(frozen=True)
class EikonalMatrixOutput:
    source: EikonalMatrixFile
    mixing: np.ndarray
    b_node: np.ndarray
    q2_node: np.ndarray
    chi_even_b: np.ndarray
    chi_odd_b: np.ndarray
    amplitude_b: tuple[np.ndarray, ...]
    amplitude_q: tuple[np.ndarray, ...]
    amplitude_zero: tuple[np.ndarray, ...]

    # Compute one diagonal pair basis opacity including even and odd exchanges
    def eigen_opacity(self, channel: tuple[int, int]) -> np.ndarray:
        pair = _pair_index(channel, self.source.nchannels)
        return self.chi_even_b[:, pair, pair] + self.chi_odd_b[:, pair, pair]

    # Project one helicity bank between physical Good Walker states
    def project_bank(self, bank: np.ndarray, final: tuple[int, int]) -> np.ndarray:
        nchannels = self.source.nchannels
        _pair_index(final, nchannels)
        initial_state = np.kron(self.mixing[0], self.mixing[0])
        final_state = np.kron(self.mixing[final[0]], self.mixing[final[1]])
        return np.einsum("i,nij,j->n", final_state, bank, initial_state, optimize=True)

    # Compute the physical spin averaged impact parameter amplitude
    def physical_impact_amplitude(self, final: tuple[int, int]) -> np.ndarray:
        return 0.25 * sum(
            (self.project_bank(self.amplitude_b[index], final) for index in _SPIN_DIAGONAL),
            start=np.zeros(self.b_node.size, dtype=complex),
        )

    # Compute the physical spin averaged momentum space amplitude
    def physical_momentum_amplitude(self, final: tuple[int, int]) -> np.ndarray:
        if self.q2_node.size == 0:
            raise ValueError(f"Eikonal matrix cache has no momentum amplitudes: {self.source.path}")
        return 0.25 * sum(
            (self.project_bank(self.amplitude_q[index], final) for index in _SPIN_DIAGONAL),
            start=np.zeros(self.q2_node.size, dtype=complex),
        )

    # Compute the exact forward physical spin averaged amplitude
    def physical_forward_amplitude(self, final: tuple[int, int]) -> complex:
        if not self.amplitude_zero:
            raise ValueError(f"Eikonal matrix cache has no forward amplitudes: {self.source.path}")
        diagonal = sum(
            (self.project_bank(self.amplitude_zero[index][None, ...], final)[0] for index in _SPIN_DIAGONAL), start=0.0j
        )
        return complex(0.25 * diagonal)


# Compute one validated pair basis index
def _pair_index(channel: tuple[int, int], nchannels: int) -> int:
    first, second = channel
    if first < 0 or second < 0 or first >= nchannels or second >= nchannels:
        raise ValueError(f"Channel {channel} is outside the N={nchannels} Good Walker basis")
    return first * nchannels + second


# Decode one complex D by D matrix bank
def _decode_matrix_bank(payload: object, field: str, points: int, dimension: int) -> np.ndarray:
    values = decode_doubles(payload, field)
    expected = 2 * points * dimension * dimension
    if values.size != expected:
        raise ValueError(f'Eikonal matrix field "{field}" has {values.size} values, expected {expected}')
    components = values.reshape(points, dimension, dimension, 2)
    return components[..., 0] + 1j * components[..., 1]


# Extract Mandelstam s from the exact cache identity
def _sqrts_from_key(key: object, path: Path) -> float:
    if not isinstance(key, str) or ";" not in key:
        raise ValueError(f"Eikonal matrix cache has no valid identity: {path}")
    try:
        mandelstam_s = float.fromhex(key.split(";", 1)[0])
    except ValueError as exc:
        raise ValueError(f"Eikonal matrix cache has invalid Mandelstam s: {path}") from exc
    if not np.isfinite(mandelstam_s) or mandelstam_s <= 0.0:
        raise ValueError(f"Eikonal matrix cache has nonpositive Mandelstam s: {path}")
    return float(np.sqrt(mandelstam_s))


# Validate filename and lightweight cache metadata
def inspect_matrix_output(path: Path) -> EikonalMatrixFile:
    match = _FILENAME.match(path.name)
    if match is None:
        raise ValueError(f"Not an eikonal matrix output filename: {path}")
    document = load_json_file(path)
    nchannels = int(match.group("n"))
    if document.get("version") != _CACHE_VERSION:
        raise ValueError(f"Unsupported eikonal matrix cache version in {path}")
    if document.get("channels") != nchannels or document.get("dimension") != nchannels**2:
        raise ValueError(f"Eikonal matrix cache dimensions disagree with its filename: {path}")
    return EikonalMatrixFile(
        path=path,
        nchannels=nchannels,
        beam1=int(match.group("b1")),
        beam2=int(match.group("b2")),
        sqrts=_sqrts_from_key(document.get("key"), path),
    )


# Discover current matrix eikonal outputs while ignoring incompatible caches
def discover_matrix_outputs(directory: Path) -> list[EikonalMatrixFile]:
    outputs: list[EikonalMatrixFile] = []
    for path in sorted(directory.glob("MATRIX_N*.json")):
        try:
            outputs.append(inspect_matrix_output(path))
        except (OSError, ValueError, json.JSONDecodeError):
            continue
    return outputs


# Read and validate the coupled channel fields exposed by this analysis
def read_matrix_output(source: EikonalMatrixFile) -> EikonalMatrixOutput:
    document = load_json_file(source.path)
    if document.get("version") != _CACHE_VERSION:
        raise ValueError(f"Unsupported eikonal matrix cache version in {source.path}")
    nchannels = source.nchannels
    dimension = nchannels**2
    if document.get("channels") != nchannels or document.get("dimension") != dimension:
        raise ValueError(f"Eikonal matrix cache dimensions changed during loading: {source.path}")

    mixing_values = decode_doubles(document.get("U"), "U")
    if mixing_values.size != nchannels**2:
        raise ValueError(f"Eikonal matrix mixing matrix has the wrong size: {source.path}")
    mixing = mixing_values.reshape(nchannels, nchannels)
    identity = mixing @ mixing.T
    if not np.allclose(identity, np.eye(nchannels), rtol=0.0, atol=5e-13):
        raise ValueError(f"Eikonal matrix Good Walker rotation is not orthogonal: {source.path}")

    b_node = decode_doubles(document.get("b_node"), "b_node")
    if b_node.size < 2 or np.any(np.diff(b_node) <= 0.0):
        raise ValueError(f"Eikonal matrix impact parameter grid is not increasing: {source.path}")
    chi_even = _decode_matrix_bank(document.get("chi_even_b"), "chi_even_b", b_node.size, dimension)
    chi_odd = _decode_matrix_bank(document.get("chi_odd_b"), "chi_odd_b", b_node.size, dimension)

    amplitude_b_data = document.get("amplitude_spin_b")
    if not isinstance(amplitude_b_data, list) or len(amplitude_b_data) != _SPIN_ENTRIES:
        raise ValueError(f"Eikonal matrix impact helicity bank is incomplete: {source.path}")
    amplitude_b = tuple(
        _decode_matrix_bank(payload, f"amplitude_spin_b[{index}]", b_node.size, dimension)
        for index, payload in enumerate(amplitude_b_data)
    )

    if bool(document.get("has_amplitude")):
        q2_node = decode_doubles(document.get("q2_node"), "q2_node")
        if q2_node.size < 2 or np.any(np.diff(q2_node) <= 0.0):
            raise ValueError(f"Eikonal matrix momentum grid is not increasing: {source.path}")
        amplitude_q_data = document.get("amplitude_spin_q")
        if not isinstance(amplitude_q_data, list) or len(amplitude_q_data) != _SPIN_ENTRIES:
            raise ValueError(f"Eikonal matrix momentum helicity bank is incomplete: {source.path}")
        amplitude_q = tuple(
            _decode_matrix_bank(payload, f"amplitude_spin_q[{index}]", q2_node.size, dimension)
            for index, payload in enumerate(amplitude_q_data)
        )
        amplitude_zero_data = document.get("amplitude_spin_zero")
        if not isinstance(amplitude_zero_data, list) or len(amplitude_zero_data) != _SPIN_ENTRIES:
            raise ValueError(f"Eikonal matrix forward helicity bank is incomplete: {source.path}")
        amplitude_zero = tuple(
            _decode_matrix_bank(payload, f"amplitude_spin_zero[{index}]", 1, dimension)[0]
            for index, payload in enumerate(amplitude_zero_data)
        )
    else:
        q2_node = np.empty(0, dtype=float)
        amplitude_q = tuple(np.empty((0, dimension, dimension), dtype=complex) for _ in range(_SPIN_ENTRIES))
        amplitude_zero = ()

    return EikonalMatrixOutput(
        source=source,
        mixing=mixing,
        b_node=b_node,
        q2_node=q2_node,
        chi_even_b=chi_even,
        chi_odd_b=chi_odd,
        amplitude_b=amplitude_b,
        amplitude_q=amplitude_q,
        amplitude_zero=amplitude_zero,
    )


# Select the newest output for one physical configuration
def select_matrix_output(
    files: list[EikonalMatrixFile], nchannels: int, beam: tuple[int, int], sqrts: float
) -> EikonalMatrixFile | None:
    matches = [
        item
        for item in files
        if item.nchannels == nchannels and (item.beam1, item.beam2) == beam and np.isclose(item.sqrts, sqrts)
    ]
    if not matches:
        return None
    return max(matches, key=lambda item: item.path.stat().st_mtime)
