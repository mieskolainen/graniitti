#!/bin/bash
#
# HepMC 3.x compilation and install script
#
# Ubuntu requirements:
# sudo apt-get install cmake
# 
# Run with:
# INSTALLPATH=$HOME/local
# source installHepMC3.sh
#
#git clone https://gitlab.cern.ch/hepmc/HepMC3.git

tar -xf HepMC3-150119.tar.gz
cd HepMC3
cmake -D HEPMC_ENABLE_ROOTIO=OFF .
cmake --build .
cmake . -DCMAKE_INSTALL_PREFIX=${INSTALLPATH}/HEPMC3
make
make install
cd ..
wait 2
rm HepMC3 -f -r
