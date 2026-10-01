#!/bin/bash
#
# HepMC 3.x compilation and install script
#
# See: https://gitlab.cern.ch/hepmc/HepMC3
#
# Run first: source setenv.sh
# 
# mikael.mieskolainen@cern.ch (2026)

set -euo pipefail

#git clone https://gitlab.cern.ch/hepmc/HepMC3.git
tar -xf HepMC3-3.3.1.tar.gz
cd HepMC3-3.3.1

# Match the upstream v0 heavy-ion serializer to its reader
# https://gitlab.cern.ch/hepmc/HepMC3/-/blob/master/src/GenHeavyIon.cc
patch -p1 < ../hepmc3-heavy-ion.patch

mkdir hepmc3-build
cd hepmc3-build

# Without dot
PYTHON_NODOT=${PYTHON_VERSION//./}

# Remove old
rm ${GRANIITTI_IO_PATH}/HEPMC3 -f -r

cmake -DCMAKE_INSTALL_PREFIX=${GRANIITTI_IO_PATH}/HEPMC3 \
	  -DHEPMC3_ENABLE_ROOTIO:BOOL=OFF \
	  -DHEPMC3_ENABLE_PROTOBUFIO:BOOL=OFF \
	  -DHEPMC3_ENABLE_TEST:BOOL=OFF \
	  -DHEPMC3_BUILD_EXAMPLES=OFF \
	  -DHEPMC3_INSTALL_INTERFACES:BOOL=ON \
	  -DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF \
	  -DHEPMC3_BUILD_DOCS:BOOL=OFF \
	  -DHEPMC3_ENABLE_PYTHON:BOOL=ON \
	  -DHEPMC3_PYTHON_VERSIONS=${PYTHON_VERSION} \
	  -DHEPMC3_Python_SITEARCH${PYTHON_NODOT}=${GRANIITTI_IO_PATH}/HEPMC3/lib/python${PYTHON_VERSION}/site-packages \
	  ../

make -j4
make install

cd ../..
sleep 3
rm HepMC3-3.3.1 -f -r
