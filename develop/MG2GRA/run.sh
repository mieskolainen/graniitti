#!/bin/sh
#
# Run the MG5 regeneration and conversion workflow through regenerate.py
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PYTHONHASHSEED=0 exec python -B "$script_dir/regenerate.py" "$@"
