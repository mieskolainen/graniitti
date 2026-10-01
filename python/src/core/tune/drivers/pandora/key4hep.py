# Prepare the selected Key4hep installation for Pandora reconstruction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re
import shlex
import subprocess
from pathlib import Path

from core.tune.drivers.pandora.runtime import command_environment, file_identity

MODULES = ("PandoraSDK", "k4geo", "k4RecTracker", "LCContent", "DDMarlinPandora")


# Read the pinned Spack environment and preserve setup diagnostics
def key4hep_environment(setup, release):
    if not setup or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", release):
        raise ValueError("Pandora requires PANDORA_KEY4HEP_SETUP and a dated PANDORA_KEY4HEP_RELEASE")
    result = subprocess.run(
        ["/bin/bash", "--noprofile", "--norc", "-c",
         'source "$1" --spack -r "$2" >&2 && /usr/bin/env -0', "pandora", str(setup), release],
        env=command_environment(), capture_output=True, text=True, timeout=60, check=False,
    )
    if result.returncode:
        raise ValueError(f"Key4hep {release} setup failed (exit {result.returncode}):\n{result.stderr.strip()}")
    return dict(entry.split("=", 1) for entry in result.stdout.split("\0") if "=" in entry)


# Extract exact package prefixes including platform, version and Spack hash
def package_prefixes(text):
    part = r"[^/\s:;'\"=\[\]]+"
    return set(re.findall(
        rf"/[^\s:;'\"=\[\]]+/key4hep/releases/\d{{4}}-\d{{2}}-\d{{2}}/{part}/{part}/{part}", text
    ))


# Allow upstream dependencies only when the selected Key4hep stack uses those exact packages
def validate_build(root, environment, release):
    root = Path(root)
    selected = environment.get("KEY4HEP_STACK", "")
    if not re.search(rf"/key4hep/releases/{re.escape(release)}/[^/]+/key4hep-stack/", selected):
        raise ValueError(f"Key4hep setup selected {selected!r}, expected release {release}")
    allowed = package_prefixes("\n".join(environment.values()))
    mismatches = []
    for module in MODULES:
        cache = root / module / "build/CMakeCache.txt"
        used = package_prefixes("\n".join(line for line in cache.read_text().splitlines() if "=" in line))
        if not used:
            raise ValueError(f"No Key4hep package paths in {cache}")
        mismatches.extend(f"{module}: {prefix}" for prefix in sorted(used - allowed))
        if not (root / module / "install").is_dir():
            raise FileNotFoundError(f"Pandora installation missing: {module}")
    if mismatches:
        raise ValueError(
            f"Pandora build packages disagree with Key4hep {release}:\n" + "\n".join(mismatches)
        )


# Generate a setup using only the pinned stack and the frozen reconstruction installs
def setup_text(setup, release):
    lines = [
        f"source {shlex.quote(setup)} --spack -r {shlex.quote(release)} || return $?",
        'pandora_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
    ]
    for module in MODULES:
        lines.extend([f'cd "$pandora_root/{module}" || return $?', "k4_local_repo || return $?"])
    lines += ['export PANDORASDK_DIR="$pandora_root/PandoraSDK"', 'export K4GEO="$pandora_root/k4geo"',
              'cd "$pandora_root" || return $?']
    return "\n".join(lines) + "\n"


# Load the same Gaudi C++ library used by k4run before a worker joins Ray
def worker_check(setup, release):
    return dict(
        command=["/bin/bash", "--noprofile", "--norc", "-c",
                 'source "$1" --spack -r "$2" && exec python -c "$3"', "pandora", str(setup), release,
                 'import ctypes, Gaudi.Main; ctypes.CDLL("libGaudiKernel.so", mode=ctypes.RTLD_GLOBAL)'],
        env=command_environment(inherit=False), timeout=60,
    )


# Supply an explicit runtime description before freezing the reconstruction inputs
def prepare(*, tunesetup, environment):
    if all("runtime" in card for card in tunesetup.datacards):
        return
    root = tunesetup.datacards[0]["pandora_dir"]
    setup, release = environment.get("PANDORA_KEY4HEP_SETUP", ""), environment.get("PANDORA_KEY4HEP_RELEASE", "")
    validate_build(root, key4hep_environment(setup, release), release)
    tunesetup.aux_param_space["worker_check"] = worker_check(setup, release)
    tunesetup.aux_param_space["pandora_runtime"] = {
        "setup": setup_text(setup, release), "dependencies": {setup: file_identity(setup)}}
