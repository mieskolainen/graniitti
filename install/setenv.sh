#!/bin/bash
#
# Set environment variables for HepMC3 and LHAPDF6
#
# For automatically loading this, add to end of your your ~/.bashrc
#
# run with: source setenv.sh

# ***
INSTALLPATH=$HOME/local
# ***

# lib64 needed on some systems
export HEPMC3SYS=${INSTALLPATH}/HEPMC3
export PATH=${HEPMC3SYS}/bin:${PATH}
export LD_LIBRARY_PATH=${HEPMC3SYS}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${HEPMC3SYS}/lib64:${LD_LIBRARY_PATH}

# lib64 needed on some systems
export LHAPDFSYS=${INSTALLPATH}/LHAPDF
export PATH=${LHAPDFSYS}/bin:${PATH}
export LD_LIBRARY_PATH=${LHAPDFSYS}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${LHAPDFSYS}/lib64:${LD_LIBRARY_PATH}

echo 'New environment variables:'
echo ''
echo ' HEPMC3SYS='$HEPMC3SYS
echo ' LHAPDFSYS='$LHAPDFSYS
echo ''
echo ' PATH='$PATH
echo ''
echo ' LD_LIBRARY_PATH='$LD_LIBRARY_PATH
