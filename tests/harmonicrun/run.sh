#!/bin/sh
#
# Run with: source ./tests/harmonic/run.sh

read -p "harmonicrun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
./bin/gr -i ./tests/harmonicrun/SH_2pi_REF.json -n 1000000
./bin/gr -i ./tests/harmonicrun/SH_2pi.json     -n 1000000

fi

# ***********************************************************************
# Fiducial cuts <ETAMIN,ETAMAX,PTMIN,PTMAX>
FIDCUTS=-0.9,0.9,0.1,100.0
# ***********************************************************************

# Maximum angular order <positive integer>
LMAX=2

# Lorentz frame (HE,CS,GJ,PG,RF)
FRAME=CS

# Mass binning <bins,min,max>
MBINS=40,0.25,2.5

# Analyze
./bin/fitharmonic -r SH_2pi_REF -i SH_2pi -c $FIDCUTS -l $LMAX -f $FRAME -M $MBINS
