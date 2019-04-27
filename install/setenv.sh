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
export PATH=${PATH}:${HEPMC3SYS}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HEPMC3SYS}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HEPMC3SYS}/lib64

# lib64 needed on some systems
export LHAPDFSYS=${INSTALLPATH}/LHAPDF
export PATH=${PATH}:${LHAPDFSYS}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LHAPDFSYS}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LHAPDFSYS}/lib64

echo 'New environment variables:'
echo ''
echo ' HEPMC3SYS='$HEPMC3SYS
echo ' LHAPDFSYS='$LHAPDFSYS
echo ''
echo ' PATH='$PATH
echo ''
echo ' LD_LIBRARY_PATH='$LD_LIBRARY_PATH
