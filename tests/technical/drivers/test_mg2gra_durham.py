# Regression tests for generated Durham MadGraph amplitudes
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
MG2GRA = ROOT / "develop" / "MG2GRA"
CARDS_ROOT = ROOT / "MG5cards"
sys.path.insert(0, str(MG2GRA))

from modules import output_layout  # noqa: E402

DURHAM_PROCESSES = (
    ("gg_ddbar", 2, 16),
    ("gg_ssbar", 2, 16),
    ("gg_ccbar", 2, 16),
    ("gg_bbbar", 2, 16),
    ("gg_uubar", 2, 16),
    ("gg_uubarg", 6, 32),
    ("gg_ddbarg", 6, 32),
    ("gg_ssbarg", 6, 32),
    ("gg_ccbarg", 6, 32),
    ("gg_bbbarg", 6, 32),
    ("gg_gg", 6, 16),
    ("gg_ggg", 24, 32),
    ("gg_gggg", 120, 64),
)


# Verify every Durham output uses the generic projected-color interface
@pytest.mark.parametrize(("process", "ncolor", "ncomb"), DURHAM_PROCESSES)
def test_durham_wrapper_conversion_dimensions(process, ncolor, ncomb):
    directory = output_layout.projection_cpp_directory("durham")
    header = (
        ROOT / output_layout.header_path(directory, f"AMP_MG5_{process}.h")
    ).read_text()
    source = (
        ROOT / output_layout.cpp_source_path(directory, f"AMP_MG5_{process}.cc")
    ).read_text()
    process_data = json.loads(
        (CARDS_ROOT / output_layout.color_path("durham", process)).read_text()
    )

    assert process_data["ncolor"] == ncolor
    assert f"const int ncomb = {ncomb};" in source
    assert f"static const int ncolor     = {ncolor};" in header
    assert f"static constexpr int nhelicity = {ncomb};" in header
    assert "#include <span>" in source
    assert "CalcColorFlowHelicity" in header
    assert "CalcColorProjectedHelicity" in header
    assert "CalcColorProjectedHelicity" in source
    assert "std::span<const std::complex<double>> jamp_view" in source
    assert "gra::BilinearProduct(projector_view, jamp_view)" in source
    assert "value += color_projectors" not in source
    assert "calculate_color_flows" in source
    assert "Contract the JAMP color basis with the full finite-Nc metric" in source
    assert "Lifson and Mattelaer" in source
    assert "mME.clear();" in source
    assert "int igood[ncomb + 1] = {};" in source
    assert "final_buffer" in header
    assert "CalcSingletHelicity" not in header
    assert "CalcSingletHelicity" not in source
