# Tests for active and fallback conda environment launchers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from submit import CAMPAIGN_DIR, campaign, schedulers
from submit.__main__ import lxplus_conda_environment
from tests.technical.support import environment as support

ROOT = Path(__file__).resolve().parents[3]
SHELL_HELPER = ROOT / "tests" / "environment.sh"


# Activate real Conda when batch jobs cannot discover the submitting user's named environment
@pytest.mark.parametrize("phase", ["init", "steer"])
@pytest.mark.parametrize("selection", ["name", "prefix"])
def test_lxplus_job_conda_activation(tmp_path, monkeypatch, phase, selection):
    envs = tmp_path / "envs"
    envs.mkdir()
    prefix = envs / "batch-graniitti"
    prefix.symlink_to(sys.prefix, target_is_directory=True)
    monkeypatch.setenv("CONDA_ENVS_PATH", str(envs))
    catalog = campaign.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml")
    if selection == "name":
        catalog["defaults"]["runtime"]["conda_env"] = prefix.name
    else:
        catalog["campaigns"]["tune-gpom-star-cms-ampfit"].setdefault("runtime", {})["conda_env"] = str(prefix)
    environment = campaign.resolve(catalog, campaign_name="tune-gpom-star-cms-ampfit")
    environment["RAY_REPO_DIR"] = str(ROOT)
    environment.update(lxplus_conda_environment(environment))
    if phase == "steer":
        rendered = schedulers.render_lxplus_job(environment, log_dir=tmp_path)
    else:
        rendered = getattr(schedulers, f"render_ray_{phase}_job")(
            environment=environment, runtime_archive=tmp_path / "runtime.tar.zst", runtime_sha256="a" * 64,
            log_dir=tmp_path, bootstrap_cache=tmp_path, coord_dir=tmp_path,
            **{f"{phase}_script": ROOT / f"submit/shell/steer_ray_{phase}_lxplus.sh"})
    assignments = next(line.split('"', 2)[1] for line in rendered.splitlines()
                       if line.startswith("environment "))
    batch_home = tmp_path / "batch-home"
    batch_home.mkdir()
    batch_environment = {"HOME": str(batch_home), "PATH": "/usr/bin:/bin", "TMPDIR": str(tmp_path)}
    batch_environment.update(dict(item.split("=", 1) for item in shlex.split(assignments)))
    missing = subprocess.run(
        [batch_environment["CONDA_EXE"], "shell.posix+json", "activate", prefix.name],
        env=batch_environment, capture_output=True, text=True)
    assert missing.returncode != 0
    assert "EnvironmentNameNotFound" in missing.stderr
    result = subprocess.run(
        ["bash", "-c", "source install/setconda_lxplus.sh || exit $?\n"
         'python -c \'import os, sys, ray; assert sys.prefix == os.environ["ICETUNE_CONDA_PREFIX"]\''],
        cwd=ROOT, env=batch_environment, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


# Reject an unavailable environment during submission before requesting batch resources
def test_lxplus_missing_conda_environment(tmp_path):
    with pytest.raises(ValueError, match="Set runtime.conda_env"):
        lxplus_conda_environment({})
    with pytest.raises(ValueError, match="Cannot resolve Conda environment"):
        lxplus_conda_environment({"ICETUNE_CONDA_ENV": str(tmp_path / "missing")})


# Do not mistake another environment with the same basename for the requested prefix
def test_shell_launcher_requires_exact_prefix(tmp_path):
    result = subprocess.run(
        ["bash", "-c", 'source "$1"; conda_environment_is_active "$2"',
         "bash", str(SHELL_HELPER), str(tmp_path / Path(sys.prefix).name)],
        cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 1, result.stdout + result.stderr


# Activate the real installation from PATH or its explicit executable
@pytest.mark.parametrize("discovery", ["path", "executable", "missing"])
def test_lxplus_conda_discovery(discovery, tmp_path):
    executable = os.environ.get("CONDA_EXE") or shutil.which("conda")
    assert executable, "Run launcher tests in the graniitti Conda environment"
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("CONDA_", "BASH_FUNC_", "ICETUNE_CONDA_"))
                   and key not in {"GRANIITTI_ENV", "BASH_ENV", "ENV"}}
    environment["PATH"] = "/usr/bin:/bin"
    environment["ICETUNE_CONDA_ENV"] = sys.prefix
    if discovery == "path":
        environment["PATH"] = f"{Path(executable).parent}:{environment['PATH']}"
    elif discovery == "executable":
        environment["CONDA_EXE"] = executable
    else:
        environment["PATH"] = str(tmp_path)
    result = subprocess.run(
        [shutil.which("bash"), "-c", "source install/setconda_lxplus.sh || exit $?\n"
         'python -c \'import os, sys; assert sys.prefix == os.environ["ICETUNE_CONDA_ENV"]\''],
        cwd=ROOT, env=environment, capture_output=True, text=True)
    if discovery == "missing":
        assert result.returncode == 64
        assert "set CONDA_EXE" in result.stderr
    else:
        assert result.returncode == 0, result.stdout + result.stderr


# Verify active graniitti tests do not launch a fragile nested conda process
def test_project_active_env(monkeypatch):
    monkeypatch.setenv("CONDA_DEFAULT_ENV", "/afs/cern.ch/user/e/example/.conda/envs/graniitti")
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    monkeypatch.setattr(support.sys, "prefix", "/usr")
    command = support.project_environment_command("./bin/gr --help")

    assert command[:2] == ["bash", "-c"]
    assert "conda" not in command


# Verify tests outside graniitti still resolve the named conda environment
def test_project_inactive_env(monkeypatch):
    monkeypatch.setenv("CONDA_DEFAULT_ENV", "base")
    monkeypatch.setenv("CONDA_PREFIX", "/opt/conda")
    monkeypatch.setattr(support.sys, "prefix", "/opt/conda")
    command = support.project_environment_command("./bin/gr --help")

    assert command[:5] == [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "graniitti",
    ]


# Verify shell launchers recognize an active absolute-prefix environment
def test_shell_active_env_prefix():
    environment = os.environ.copy()
    environment["CONDA_DEFAULT_ENV"] = "/afs/cern.ch/user/e/example/.conda/envs/graniitti"
    environment.pop("CONDA_PREFIX", None)
    command = "\n".join(
        [
            "conda() { return 91; }",
            f"source {shlex.quote(str(SHELL_HELPER))}",
            "activate_conda_environment graniitti",
            "run_in_graniitti /bin/true",
        ]
    )

    result = subprocess.run(
        ["bash", "-c", command],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
