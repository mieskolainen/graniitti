# Validate the chi_cJ radiative multipole conversion tool
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import cmath
import json
import math
import pathlib
import subprocess
import sys

import pytest
import sympy as sp
from sympy.physics.wigner import clebsch_gordan

from develop.tools import jpsi_gamma_helicity

ROOT = pathlib.Path(__file__).resolve().parents[3]


# Check pure E1 gives the standard orthogonal multipole-helicity coefficients
@pytest.mark.parametrize(
    ("spin", "expected"),
    [
        (0, (1.0,)),
        (1, (1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0))),
        (
            2,
            (
                math.sqrt(1.0 / 10.0),
                math.sqrt(3.0 / 10.0),
                math.sqrt(3.0 / 5.0),
            ),
        ),
    ],
)
def test_pure_e1_helicity_coefficients(spin, expected):
    result = jpsi_gamma_helicity.convert(spin)
    assert result.helicities == pytest.approx(expected)
    assert sum(value * value for value in result.helicities) == pytest.approx(1.0)


# Derive every multipole coefficient and relative helicity phase from angular momentum addition
# [REFERENCE: Artuso et al., arXiv:0910.0046, Eq. (4)]
@pytest.mark.parametrize("spin,m2,e3", [(0, 0.0, 0.0), (1, 0.8, 0.0), (1, -0.3, 0.0),
                                     (2, -0.8, 0.5), (2, 0.2, -0.7)])
def test_multipoles_complex_rows_clebsch_gordan(spin, m2, e3):
    multipoles = (math.sqrt(1.0 - m2**2 - e3**2), m2, e3)[:spin + 1]
    matrix = sp.Matrix([[sp.sqrt(sp.Rational(2 * order + 1, 2 * spin + 1))
                        * clebsch_gordan(order, 1, spin, 1, nu - 1, nu)
                        for order in range(1, spin + 2)] for nu in range(spin + 1)])
    assert sp.simplify(matrix.T * matrix) == sp.eye(spin + 1)
    expected = [float(value) for value in matrix * sp.Matrix(multipoles)]
    result = jpsi_gamma_helicity.convert(spin, m2, e3)
    assert result.helicities == pytest.approx(expected, abs=1e-14)
    reference = expected[2 if spin == 2 else 0]
    for first, second, magnitude, phase in result.rows:
        assert abs(second) == 1
        assert magnitude * cmath.exp(1j * phase) == pytest.approx(expected[abs(first - second)] / reference, abs=1e-14)


# Check the CLI preserves the sign of a helicity amplitude for an M2 admixture
def test_converter_command_line_json():
    output = subprocess.check_output(
        [
            sys.executable,
            str(ROOT / "develop" / "tools" / "jpsi_gamma_helicity.py"),
            "--spin",
            "1",
            "--m2",
            "0.8",
            "--json",
        ],
        text=True,
    )
    payload = json.loads(output)
    assert payload[0]["spin"] == 1
    helicities = list(payload[0]["helicities"].values())
    assert helicities == pytest.approx([1.4 / math.sqrt(2), -0.2 / math.sqrt(2)])
    rows = payload[0]["helicity"]
    assert rows[1][2] * cmath.exp(1j * rows[1][3]) == pytest.approx(-1.0 / 7.0)
