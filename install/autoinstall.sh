#!/bin/bash
#
# Automatic installation of HEPMC3 and LHAPDF6 libraries
#
# Ubuntu requirements:
#   sudo apt install cmake g++ python3-dev
#
# Run with: bash autoinstall.sh under the install folder
# 
# mikael.mieskolainen@cern.ch (2026)

set -euo pipefail

# Check that the selected C++ compiler produces runnable programs
check_cxx() {
    local cxx=${CXX:-c++}
    local probe

    if ! command -v "${cxx}" >/dev/null 2>&1; then
        echo "** C++ compiler not found: ${cxx} **" >&2
        return 1
    fi

    probe=$(mktemp "${TMPDIR:-/tmp}/graniitti-cxx.XXXXXX")
    if ! printf 'int main() { return 0; }\n' | "${cxx}" -x c++ - -o "${probe}"; then
        rm -f "${probe}"
        echo "** C++ compiler test build failed: ${cxx} **" >&2
        return 1
    fi
    if ! "${probe}"; then
        rm -f "${probe}"
        echo "** C++ compiler test program cannot run: ${cxx} **" >&2
        echo '** Update the graniitti Conda environment from environment.yml **' >&2
        return 1
    fi
    rm -f "${probe}"
}

INSTALL_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "${INSTALL_DIR}"

if [[ -n "${GRANIITTI_ENV:-}" ]]; then
    echo '** GRANIITTI environment is set **'
    check_cxx

    read -p "Do you want to install HepMC3 and LHAPDF6 to ${GRANIITTI_IO_PATH} (old installation will be removed)? [y/n]" -n 1 -r

    echo # New line
    if [[ ${REPLY} =~ ^[Yy]$ ]]; then

        # Install
        bash installHEPMC3.sh
        bash installLHAPDF6.sh
        bash installPDFSET.sh

        echo # New line
        echo "Installation of libraries is done."
    fi
else
    echo '** GRANIITTI environment setup not found **'
    echo # new line
    echo 'Execute: source install/setenv.sh before running this command'
fi
