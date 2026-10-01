# Shared conda environment selection for subprocess-based tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os
import pathlib
import sys


# Compute whether pytest already runs inside the requested conda environment
def conda_environment_is_active(environment_name):
    environment_values = {
        os.environ.get("CONDA_DEFAULT_ENV", ""),
        os.environ.get("CONDA_PREFIX", ""),
        sys.prefix,
    }
    environment_names = {
        pathlib.Path(value.rstrip("/")).name for value in environment_values if value
    }
    return environment_name in environment_values or environment_name in environment_names


# Compute whether pytest already runs inside the graniitti conda environment
def graniitti_environment_is_active():
    return conda_environment_is_active("graniitti")


# Build a project shell command without nesting an already-active conda environment
def project_environment_command(command):
    shell_command = [
        "bash",
        "-c",
        f"source install/setenv.sh && {command}",
    ]
    if graniitti_environment_is_active():
        return shell_command
    return [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "graniitti",
        *shell_command,
    ]
