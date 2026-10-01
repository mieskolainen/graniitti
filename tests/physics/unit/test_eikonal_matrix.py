# Test coupled channel eikonal matrix output analysis
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import base64
import json
import zlib
from pathlib import Path

import numpy as np
import pytest
from core.tune.drivers.graniitti.eikonal import (
    discover_matrix_outputs,
    inspect_matrix_output,
    read_matrix_output,
    select_matrix_output,
)

from tests.physics.studies.elastic.analysis import analyze_impact, analyze_momentum
from tests.physics.studies.elastic.analysis.steering import scan_outputs


# Encode doubles using the production cache byte convention
def _encode(values: np.ndarray) -> str:
    raw = np.asarray(values, dtype=">f8").tobytes()
    return base64.b64encode(zlib.compress(raw, level=9)).decode("ascii")


# Encode one complex matrix bank in row major order
def _encode_bank(bank: np.ndarray) -> str:
    components = np.stack((bank.real, bank.imag), axis=-1)
    return _encode(components.reshape(-1))


# Construct one complete N=2 cache with analytically projectable operators
def _write_cache(path: Path) -> None:
    nchannels = 2
    dimension = nchannels**2
    theta = 0.37
    mixing = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    b_node = np.array([0.0, 0.5, 1.0])
    q2_node = np.array([0.1, 0.4])
    pair_operator = np.diag([1.0 + 0.2j, 2.0 - 0.3j, 3.0 + 0.4j, 5.0 + 0.7j])
    pair_operator[0, 3] = pair_operator[3, 0] = 0.4 - 0.3j

    chi_even = np.zeros((b_node.size, dimension, dimension), dtype=complex)
    chi_odd = np.zeros_like(chi_even)
    for point in range(b_node.size):
        chi_even[point] = np.diag(np.arange(dimension) + point + 1j * (2.0 + point))
        chi_odd[point] = np.diag(0.25 * np.arange(dimension) + 0.5j)

    amplitude_b = []
    amplitude_q = []
    amplitude_zero = []
    screening_b = []
    screening_q = []
    screening_zero = []
    crossing_even_q = []
    crossing_even_zero = []
    crossing_odd_q = []
    crossing_odd_zero = []
    for helicity in range(16):
        impact_scale = 2.0 + helicity
        momentum_scale = 10.0 + 2.0 * helicity
        amplitude_b.append(
            _encode_bank(np.stack([(impact_scale + point) * pair_operator for point in range(3)]))
        )
        amplitude_q.append(
            _encode_bank(np.stack([(momentum_scale + point) * pair_operator for point in range(2)]))
        )
        amplitude_zero.append(_encode_bank(((momentum_scale - 1.0) * pair_operator)[None, ...]))
        screening_b.append(amplitude_b[-1])
        screening_q.append(amplitude_q[-1])
        screening_zero.append(amplitude_zero[-1])
        crossing_even_q.append(amplitude_q[-1])
        crossing_even_zero.append(amplitude_zero[-1])
        crossing_odd_q.append(amplitude_q[-1])
        crossing_odd_zero.append(amplitude_zero[-1])

    document = {
        "version": 1,
        "key": f"{float(13000.0**2).hex()};2212;2212;test",
        "channels": nchannels,
        "dimension": dimension,
        "has_amplitude": True,
        "max_singular_value": 1.0,
        "U": _encode(mixing.reshape(-1)),
        "b_node": _encode(b_node),
        "q2_node": _encode(q2_node),
        "chi_even_b": _encode_bank(chi_even),
        "chi_odd_b": _encode_bank(chi_odd),
        "chi_cut_b": _encode_bank(chi_even),
        "amplitude_spin_b": amplitude_b,
        "amplitude_spin_q": amplitude_q,
        "amplitude_spin_zero": amplitude_zero,
        "crossing_even_spin_q": crossing_even_q,
        "crossing_even_spin_zero": crossing_even_zero,
        "crossing_odd_spin_q": crossing_odd_q,
        "crossing_odd_spin_zero": crossing_odd_zero,
        "screening_spin_b": screening_b,
        "screening_spin_q": screening_q,
        "screening_spin_zero": screening_zero,
    }
    path.write_text(json.dumps(document), encoding="utf-8")


# Expand the complex coupled channel projection as independent trigonometric polynomials
def _expected_projection(final: tuple[int, int]) -> complex:
    c, s = np.cos(0.37), np.sin(0.37)
    a, b, d, e, transition = 1.0 + 0.2j, 2.0 - 0.3j, 3.0 + 0.4j, 5.0 + 0.7j, 0.4 - 0.3j
    return {
        (0, 0): a * c**4 + (b + d + 2 * transition) * c**2 * s**2 + e * s**4,
        (0, 1): c * s * ((b - a + transition) * c**2 + (e - d - transition) * s**2),
        (1, 0): c * s * ((d - a + transition) * c**2 + (e - b - transition) * s**2),
        (1, 1): (a - b - d + e) * c**2 * s**2 + transition * (c**4 + s**4),
    }[final]


# Check exact decoding and physical Good Walker projections
def test_matrix_state_projection(tmp_path: Path) -> None:
    path = tmp_path / "MATRIX_N2_2212_2212_123.json"
    _write_cache(path)
    source = inspect_matrix_output(path)
    output = read_matrix_output(source)

    assert source.nchannels == 2
    assert source.sqrts == pytest.approx(13000.0)
    expected_mixing = np.array(
        [
            [np.cos(0.37), np.sin(0.37)],
            [-np.sin(0.37), np.cos(0.37)],
        ]
    )
    assert np.array_equal(output.mixing, expected_mixing)
    assert np.allclose(output.mixing @ output.mixing.T, np.eye(2), rtol=0.0, atol=2e-16)
    assert np.array_equal(output.b_node, [0.0, 0.5, 1.0])
    assert np.array_equal(output.q2_node, [0.1, 0.4])

    expected_opacity = np.array([2.5 + 2.5j, 3.5 + 3.5j, 4.5 + 4.5j])
    assert np.array_equal(output.eigen_opacity((1, 0)), expected_opacity)
    for final_state in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        projection = _expected_projection(final_state)
        assert projection != pytest.approx(0.0, abs=1e-14)
        assert np.allclose(
            output.physical_impact_amplitude(final_state),
            projection * np.array([9.5, 10.5, 11.5]),
            rtol=0.0,
            atol=2e-14,
        )
        assert np.allclose(
            output.physical_momentum_amplitude(final_state),
            projection * np.array([25.0, 26.0]),
            rtol=0.0,
            atol=5e-14,
        )
        assert output.physical_forward_amplitude(final_state) == pytest.approx(
            projection * 24.0,
            abs=5e-14,
        )


# Check discovery, model selection and incompatible cache handling
def test_matrix_current_cache(tmp_path: Path) -> None:
    current = tmp_path / "MATRIX_N2_2212_2212_123.json"
    _write_cache(current)
    incompatible = tmp_path / "MATRIX_N2_2212_2212_456.json"
    old_document = json.loads(current.read_text(encoding="utf-8"))
    old_document["version"] = 7
    incompatible.write_text(json.dumps(old_document), encoding="utf-8")

    files = discover_matrix_outputs(tmp_path)
    assert [item.path for item in files] == [current]
    assert select_matrix_output(files, 2, (2212, 2212), 13000.0) == files[0]
    assert select_matrix_output(files, 1, (2212, 2212), 13000.0) is None
    with pytest.raises(ValueError, match="Unsupported eikonal matrix cache version"):
        inspect_matrix_output(incompatible)


# Check both plotters use the exact scan output despite another matching cache
def test_elastic_scan_selects_recorded_output(tmp_path: Path) -> None:
    selected = tmp_path / "MATRIX_N2_2212_2212_123.json"
    unrelated = tmp_path / "MATRIX_N2_2212_2212_456.json"
    _write_cache(selected)
    _write_cache(unrelated)
    scan_dir = tmp_path / "pp"
    scan_dir.mkdir()
    log = scan_dir / "scan.log"
    log.write_text(f"Loaded matrix eikonal cache: {selected}\n[xscan: done]\n")
    files = scan_outputs(tmp_path, "double", "pp")
    assert [source.path for source in files] == [selected]
    momentum = analyze_momentum.load_channel_data(files, [13000.0], [1.0], "pp", 2, (0, 0))
    impact = analyze_impact.load_channel_data(files, [13000.0], "pp", 2, (0, 0))
    assert momentum[0][2][:, 1] == pytest.approx(
        _expected_projection((0, 0)).real * np.array([24.0, 25.0, 26.0])
    )
    assert impact[0][1][:, 3] == pytest.approx(
        _expected_projection((0, 0)).real * np.array([9.5, 10.5, 11.5])
    )
    with pytest.raises(ValueError, match="disagrees"):
        scan_outputs(tmp_path, "single", "pp")
    log.write_text(f"Loaded matrix eikonal cache: {selected}\n")
    with pytest.raises(ValueError, match="incomplete"):
        scan_outputs(tmp_path, "double", "pp")


# Check both plotters reject missing or ambiguous scan energies
@pytest.mark.parametrize("count", [0, 2])
def test_elastic_scan_requires_one_output(tmp_path: Path, count: int) -> None:
    path = tmp_path / "MATRIX_N2_2212_2212_123.json"
    _write_cache(path)
    files = [inspect_matrix_output(path)] * count
    with pytest.raises(ValueError, match="exactly one steered output"):
        analyze_momentum.load_channel_data(files, [13000.0], [1.0], "pp", 2, (0, 0))
    with pytest.raises(ValueError, match="exactly one steered output"):
        analyze_impact.load_channel_data(files, [13000.0], "pp", 2, (0, 0))


# Check that a nonorthogonal Good Walker rotation is rejected
def test_matrix_output_rejects_nonorthogonal_mixing(tmp_path: Path) -> None:
    path = tmp_path / "MATRIX_N2_2212_2212_123.json"
    _write_cache(path)
    document = json.loads(path.read_text(encoding="utf-8"))
    document["U"] = _encode(np.array([1.0, 0.0, 0.5, 1.0]))
    path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(ValueError, match="rotation is not orthogonal"):
        read_matrix_output(inspect_matrix_output(path))


# Integrate elastic predictions over measured bins instead of using their centers
def test_elastic_bin_average_and_range():
    from tests.physics.studies.elastic.analysis.analyze_momentum import bin_average

    x = np.array([0.0, 0.5, 1.0, 2.0])
    assert bin_average(x, 1.0 + 2.0 * x, np.array([0.25, 0.75, 1.5])) == pytest.approx([2.0, 3.25])
    with pytest.raises(ValueError, match="inside the ordered amplitude grid"):
        bin_average(x, x, np.array([0.0, 3.0]))
