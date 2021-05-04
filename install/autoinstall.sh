#!/bin/bash
#
# Automatic installation (tested on Ubuntu 20.1 LTS)
#
# Ubuntu requirements:
#   sudo apt install cmake g++ python3-dev
#
# Run with: source autoinstall.sh

# ****
PYTHON_VERSION=3.8
INSTALLPATH=$HOME/local
# ****

read -p "Do you want to install HepMC3 and LHAPDF6 to $INSTALLPATH (old installation will be removed)? [y/n]" -n 1 -r
echo # New line
if [[ $REPLY =~ ^[Yy]$ ]]
then

# Install
source installHEPMC3.sh
source installLHAPDF6.sh
source installPDFSET.sh

echo # New line
echo "Installation done. Now run: source setenv.sh"

fi
