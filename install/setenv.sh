#!/bin/bash
# 
# Set environment variables for HepMC3 and LHAPDF6, and Python tools
# 
# For automatically loading this, add to the end of your ~/.bashrc
# 
# Run with:
#   source install/setenv.sh
#   or
#   CONDA_CXX_ABI=1 source install/setenv.sh (Conda C++ runtime libraries first)
# 
# mikael.mieskolainen@cern.ch (2026)

# ***
# MODIFY THIS (HepMC3, LHAPDF6 will be installed here)
export GRANIITTI_IO_PATH="$HOME/local"
# ***

# Determine the available Python runtime version
if command -v python >/dev/null 2>&1; then
    PYTHON_VERSION=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
else
    printf '\033[31m%s\033[0m\n' \
        '** GRANIITTI: Python executable not found in PATH, modern analysis tools will not work, continuing. **'
    if command -v python3 >/dev/null 2>&1; then
        PYTHON_VERSION=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
    else
        PYTHON_VERSION=
    fi
fi

if [[ -n "${CONDA_CXX_ABI+x}" && "${CONDA_CXX_ABI}" != "0" && "${CONDA_CXX_ABI}" != "1" ]]; then
    echo "Invalid CONDA_CXX_ABI=${CONDA_CXX_ABI}; using automatic default"
    unset CONDA_CXX_ABI
fi

if [[ -z "${CONDA_PREFIX:-}" || ! -d "${CONDA_PREFIX}/lib" ]]; then
    export CONDA_CXX_ABI=0
else
    # Keep Conda libraries opt-in because they can shadow the ROOT/LCG C++ runtime
    export CONDA_CXX_ABI=${CONDA_CXX_ABI:-0}
fi

if [[ -n "${GRANIITTI_ENV:-}" ]]; then
    printf '\033[31m%s\033[0m\n' \
        '** GRANIITTI environment already set. Run with: GRANIITTI_ENV=; source install/setenv.sh to reset **'
    echo ''

else
    export PYTHON_VERSION=${PYTHON_VERSION}

    # lib64 needed on some systems, keep these as first in PATH for priority
    export HEPMC3SYS=${GRANIITTI_IO_PATH}/HEPMC3
    export PATH=${HEPMC3SYS}/bin:${PATH}
    export LD_LIBRARY_PATH=${HEPMC3SYS}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
    export LD_LIBRARY_PATH=${HEPMC3SYS}/lib64:${LD_LIBRARY_PATH}

    # lib64 needed on some systems, keep these as first in PATH for priority
    export LHAPDFSYS=${GRANIITTI_IO_PATH}/LHAPDF
    export PATH=${LHAPDFSYS}/bin:${PATH}
    export LD_LIBRARY_PATH=${LHAPDFSYS}/lib:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=${LHAPDFSYS}/lib64:${LD_LIBRARY_PATH}
    
    if [[ -n "${PYTHON_VERSION}" ]]; then
        export PYTHONPATH=${HEPMC3SYS}/lib/python${PYTHON_VERSION}/site-packages${PYTHONPATH:+:${PYTHONPATH}}
    else
        export PYTHONPATH=${PYTHONPATH:-}
    fi
    
    export GRANIITTI_ENV=True
fi

# Use Python tools directly from this checkout
_GRANIITTI_PYTHON="$(cd "$(dirname "${BASH_SOURCE[0]}")/../python" && pwd)"
if [[ -d "${_GRANIITTI_PYTHON}/src" ]]; then
    _GRANIITTI_PYTHON="${_GRANIITTI_PYTHON}/src"
fi
case ":${PYTHONPATH:-}:" in
    *":${_GRANIITTI_PYTHON}:"*) ;;
    *) export PYTHONPATH="${_GRANIITTI_PYTHON}${PYTHONPATH:+:${PYTHONPATH}}" ;;
esac
unset _GRANIITTI_PYTHON

# Optionally keep Conda runtime libraries first for packages that require them
if [[ "${CONDA_CXX_ABI}" == "0" ]]; then
    printf '\033[33m%s\033[0m\n' \
        'Skipping Conda runtime libraries priority for LD_LIBRARY_PATH [CONDA_CXX_ABI=0]'
    if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/lib" ]]; then
        _CONDALIB="${CONDA_PREFIX}/lib"
        LD_LIBRARY_PATH=":${LD_LIBRARY_PATH:-}:"
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH//:${_CONDALIB}:/:}"
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH#:}"
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH%:}"
        export LD_LIBRARY_PATH
        unset _CONDALIB
    fi
elif [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/lib" ]]; then
    echo 'Setting Conda runtime libraries priority for LD_LIBRARY_PATH [CONDA_CXX_ABI=1]'
    _CONDALIB="${CONDA_PREFIX}/lib"
    LD_LIBRARY_PATH=":${LD_LIBRARY_PATH:-}:"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH//:${_CONDALIB}:/:}"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH#:}"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH%:}"
    export LD_LIBRARY_PATH="${_CONDALIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    unset _CONDALIB
fi

echo '' 
printf '\033[32m%s\033[0m\n' '** GRANIITTI: Environment variables set **'
echo ''
echo ''
echo ' HEPMC3SYS='${HEPMC3SYS:-}
echo ' LHAPDFSYS='${LHAPDFSYS:-}
echo ' ROOTSYS='${ROOTSYS:-}
echo ' LIBTORCHSYS='${LIBTORCHSYS:-}
echo ''
echo ' GRANIITTI_IO_PATH='${GRANIITTI_IO_PATH:-}
echo ' CONDA_CXX_ABI='${CONDA_CXX_ABI:-}
echo ''
echo ' PATH='${PATH:-}
echo ''
echo ' LD_LIBRARY_PATH='${LD_LIBRARY_PATH:-}
echo ''
echo ' PYTHONPATH='${PYTHONPATH:-}
echo ' PYTHON_VERSION='${PYTHON_VERSION:-}
echo ''

printf '\033[33m%s\033[0m\n' '** System memory limits **'
echo ''
# Keep existing limits when restricted runners do not permit changes
ulimit -s unlimited 2>/dev/null || echo 'Stack limit unchanged (permission denied)'
ulimit -u 65536 2>/dev/null || echo 'Process limit unchanged (permission denied)'
ulimit -a
echo ''
