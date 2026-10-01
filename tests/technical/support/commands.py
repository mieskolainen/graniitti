# Helper functions for unit & global tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import subprocess
import sys


# Execute one shell command and return its combined output
def execute_one(command):
    result = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    sys.stdout.write(result.stdout)
    assert result.returncode == 0, (
        f"Command failed with exit code {result.returncode}: {command}\n{result.stdout}"
    )
    return result.stdout


# Execute one or more commands and require the expected completion marker
def execute(cmd, expect="[gr: done]"):
    commands = cmd if isinstance(cmd, list) else [cmd]
    outputs = []
    for command in commands:
        print(f"{__name__}: executing: {command}")
        output = execute_one(command)
        if expect is not None:
            assert expect in output, f"Command output is missing {expect!r}: {command}\n{output}"
        outputs.append(output)
    return outputs if isinstance(cmd, list) else outputs[0]
