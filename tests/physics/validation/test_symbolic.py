# Pytest entry points for exact symbolic study checks
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import runpy
from pathlib import Path

import pytest

SYMBOLIC_ROOT = Path(__file__).resolve().parents[1] / "symbolic"
SYMBOLIC_SCRIPTS = tuple(sorted(SYMBOLIC_ROOT.glob("*.py")))


# Run every exact symbolic derivation discovered in the symbolic directory
@pytest.mark.parametrize("script", SYMBOLIC_SCRIPTS, ids=lambda path: path.stem)
def test_symbolic_derivation(script):
    namespace = runpy.run_path(str(script))
    namespace["main"]()
