#!/usr/bin/env bash
# Format the requested C++ files
#
# Run with: bash develop/format.sh sourcefile.cc

set -euo pipefail

if (( $# == 0 )); then
    echo "Usage: bash develop/format.sh sourcefile.cc [...]" >&2
    exit 2
fi

FORMAT="clang-format"
if [[ -x /usr/local/bin/clang-format ]]; then
    FORMAT="/usr/local/bin/clang-format"
fi

printf 'clang-format processing %s files ...\n' "$#"
"$FORMAT" -fallback-style=none -i "$@"
echo "clang-format done!"
