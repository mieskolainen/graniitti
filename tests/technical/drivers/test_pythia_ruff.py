# Test for repository Python Ruff compliance
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


# Require the repository Python surface to pass the Ruff policy
def test_repository_python_passes_ruff():
    executable = shutil.which("ruff")
    assert executable is not None, "Ruff is required; install requirements.txt"
    result = subprocess.run(
        [executable, "check", "."],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
