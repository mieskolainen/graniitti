#!/bin/bash
#
# HepMC 3.x compilation and install script
#
# Ubuntu requirements:
#   sudo apt install cmake g++ python3-dev
#
# See: https://gitlab.cern.ch/hepmc/HepMC3
# 
# Run with:
#   
#   PYTHON_VERSION=3.8
#   INSTALLPATH=$HOME/local
#   source installHEPMC3.sh
#

#git clone https://gitlab.cern.ch/hepmc/HepMC3.git
tar -xf HepMC3-3.2.3.tar.gz
cd HepMC3-3.2.3

mkdir hepmc3-build
cd hepmc3-build

# Remove old
rm ${INSTALLPATH}/HEPMC3 -f -r

cmake -DHEPMC3_ENABLE_ROOTIO:BOOL=OFF -DHEPMC3_ENABLE_TEST:BOOL=OFF  \
-DHEPMC3_INSTALL_INTERFACES:BOOL=ON -DHEPMC3_ENABLE_PYTHON:BOOL=ON -DHEPMC3_PYTHON_VERSIONS=${PYTHON_VERSION}  \
-DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF -DHEPMC3_BUILD_DOCS:BOOL=OFF  \
-DCMAKE_INSTALL_PREFIX=${INSTALLPATH}/HEPMC3   \
-DHEPMC3_Python_SITEARCH38=${INSTALLPATH}/HEPMC3/lib/python${PYTHON_VERSION}/site-packages \
../

make -j4
make install

cd ../..
sleep 3
rm HepMC3-3.2.3 -f -r

