#!/bin/sh
#
# install mcviz HepMC2 event visualization
#
# run with: source install.sh

git clone git://github.com/mcviz/mcviz mcviz
cd mcviz
python2.7 mcviz_bootstrap.py
cd ..

