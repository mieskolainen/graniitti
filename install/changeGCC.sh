#!/bin/sh
#
# Change Ubuntu 16.04 GCC version
# Run with: source changeGCC.sh
#
# (basecode from gist.github.com/jlblancoc)
#
# Tested to change succesfully between 5 <-> 7

GCC_VER=7
PRIORITY=10

read -p "Enter GCC VERSION (e.g. 5 or 7): "  GCC_VER

sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-$GCC_VER $PRIORITY \
                         --slave /usr/bin/g++ g++ /usr/bin/g++-$GCC_VER
sudo update-alternatives --config gcc
gcc --version
g++ --version

# Map all toolchain to new version
ls -la /usr/bin/ | grep -oP "[\S]*(gcc|g\+\+)(-[a-z]+)*[\s]" | xargs sudo bash -c 'for link in ${@:1}; do ln -s -f "/usr/bin/${link}-${0}" "/usr/bin/${link}"; done' $GCC_VER

