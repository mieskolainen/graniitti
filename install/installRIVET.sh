#!/bin/bash
#
# Install the full RIVET environment
#
#
# Ubuntu requirements:
#   sudo apt install cmake g++ python3-dev
#
# First make sure HEPMC3 is installed
# 
# Run with:
#  source installRIVET.sh
# 

# For HEPMC3!
source setenv.sh


# Make temp dir
mkdir temp
cp rivet-bootstrap ./temp
cd temp

# Install
# Note we force HEPMC3 inside rivet-bootstrap (see comments within the file)
INSTALL_PREFIX=$HOME/local/rivet MAKE="make -j8" ./rivet-bootstrap

cd ..
rm temp -f -r

# Set the environment
source $HOME/local/rivet/rivetenv.sh



#  wget https://gitlab.com/hepcedar/rivetbootstrap/raw/3.1.4/rivet-bootstrap
#  chmod +x rivet-bootstrap