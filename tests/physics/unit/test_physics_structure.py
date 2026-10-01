# Check physics test selection markers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[3]
ICEPACK_ROOT = ROOT / "icepack"
VALIDATION_ROOT = ROOT / "tests" / "physics" / "validation"


# Select publication comparisons and numerical closures without classifying execution alone as agreement
@pytest.mark.parametrize(
    ("dataset", "expected"),
    [
        ("UPC/PHOTOPROD/CMS_2908607/phi_coherent", {"physics", "literature"}),
        ("MC/pipi", {"physics", "closure"}),
        ("HARDPOM/z", {"physics"}),
    ],
)
def test_physics_markers(dataset, expected):
    from tests.conftest import pytest_collection_modifyitems

    markers = []
    item = SimpleNamespace(
        path=VALIDATION_ROOT / "test_icepacks.py",
        callspec=SimpleNamespace(params={"dataset_path": ICEPACK_ROOT / dataset / "dataset.json"}),
        originalname="test_icepack_event_flow",
        add_marker=markers.append,
    )
    config = SimpleNamespace(getoption=lambda option: option == "--run-physics")
    pytest_collection_modifyitems(config, [item])
    assert {marker.name for marker in markers} == expected
