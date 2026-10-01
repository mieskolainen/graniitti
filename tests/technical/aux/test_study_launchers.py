# Check the shared study launcher
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
COMMON = "source tests/physics/studies/common.sh"


# Resolve one study action without invoking a study executable
def resolve_action(action: str | None) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    if action is None:
        environment.pop("STUDY_ACTION", None)
    else:
        environment["STUDY_ACTION"] = action
    return subprocess.run(
        [
            "bash",
            "-c",
            f'{COMMON}; study_resolve_action && printf "%s\\n" "$STUDY_ACTION"',
        ],
        cwd=ROOT,
        env=environment,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=5,
        check=False,
    )


# Accept both supported explicit actions
def test_action_accepts_both_and_analyze():
    for action in ("both", "analyze"):
        result = resolve_action(action)
        assert result.returncode == 0
        assert result.stdout == f"{action}\n"


# Default a noninteractive launcher to generation and analysis
def test_action_noninteractive_default_both():
    result = resolve_action(None)
    assert result.returncode == 0
    assert result.stdout == "both\n"


# Reject unknown actions before a study can start
def test_study_action_rejects_unknown_value():
    result = resolve_action("generate")
    assert result.returncode == 64


# Write one executable which prints each received argument on its own line
def write_argument_printer(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$@"\n', encoding="utf-8")
    path.chmod(0o755)


# Require generator wrappers to replace card defaults with shared runtime flags
def test_study_gr_runtime_arguments(tmp_path):
    write_argument_printer(tmp_path / "bin" / "gr")
    result = subprocess.run(
        [
            "bash",
            "-c",
            (
                f'{COMMON}; STUDY_REPO_ROOT="$1"; NEVENTS=23; '
                "WEIGHTED=true; LOOPSCREEN=false; "
                "study_gr -i card.json -n 5 -w false -l true"
            ),
            "study-gr-test",
            str(tmp_path),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [
        "-i",
        "card.json",
        "--NEVENTS",
        "23",
        "--WEIGHTED",
        "true",
        "--LOOPSCREEN",
        "false",
    ]


# Require both C++ analyzers to inherit the shared event record limit
@pytest.mark.parametrize(
    ("executable", "function"),
    (("analyze", "study_analyze"), ("fitharmonic", "study_fitharmonic")),
)
def test_study_analyzer_runtime_arguments(tmp_path, executable, function):
    write_argument_printer(tmp_path / "bin" / executable)
    result = subprocess.run(
        [
            "bash",
            "-c",
            (
                f'{COMMON}; STUDY_REPO_ROOT="$1"; NEVENTS=37; '
                f"{function} --card analysis.json --maximum 9"
            ),
            "study-analyzer-test",
            str(tmp_path),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [
        "--card",
        "analysis.json",
        "--maximum",
        "37",
    ]


# Require every study launcher and aggregate runner to parse as Bash
def test_all_study_launchers_pass_bash_syntax():
    study_root = ROOT / "tests" / "physics" / "studies"
    launchers = sorted(study_root.rglob("run.sh"))
    launchers.append(study_root / "run_all.sh")
    launchers.append(study_root / "run_plots.sh")
    for launcher in launchers:
        result = subprocess.run(
            ["bash", "-n", str(launcher)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, f"{launcher}: {result.stderr}"
