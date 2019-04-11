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

# Lorentz frame (HE,CS,GJ,PG,SR)
FRAME=HE

# System kinematic variables binning <bins,min,max>
MBINS=40,0.25,1.5
PBINS=1,0.0,1.0
YBINS=40,-0.9,0.9

# PARAMETERS
LMAX=4
REMOVEODD=true
REMOVENEGATIVE=true
LAMBDA=0.001
EML=false

# Fast simulation
FS=true

# Analyze
#./bin/fitharmonic -r SH_2pi_REF -i ALICE7_2pi.csv \
./bin/fitharmonic -r SH_2pi_REF -i SH_2pi \
-c $FIDCUTS -f \
$FRAME -l $LMAX -o $REMOVEODD -n $REMOVENEGATIVE -a $LAMBDA -e $EML \
-M $MBINS -P $PBINS -Y $YBINS \
-s $FS -X 100000

# Implement 2D harmonic plots (M,Pt)
# ...
