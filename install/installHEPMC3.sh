#!/bin/bash
#
# HepMC 3.x compilation and install script
#
# Ubuntu requirements:
# sudo apt-get install cmake
#
# Run with:
#
# INSTALLPATH=$HOME/local
# source installHEPMC3.sh
#
#git clone https://gitlab.cern.ch/hepmc/HepMC3.git

#tar -xf HepMC3-130319.tar.gz
tar -xf HepMC3-3.2.0.tar.gz

cd HepMC3-3.2.0

cmake -DHEPMC3_ENABLE_SEARCH=ON -DHEPMC3_ENABLE_PYTHON=OFF -DHEPMC3_ENABLE_ROOTIO=OFF .
cmake --build .
cmake . -DCMAKE_INSTALL_PREFIX=${INSTALLPATH}/HEPMC3
make -j4
make test
make install
cd ..
wait 2
rm HepMC3-3.2.0 -f -r

