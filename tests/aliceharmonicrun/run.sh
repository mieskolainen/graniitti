#!/bin/sh
#
# Run with: source ./tests/harmonic/run.sh

read -p "aliceharmonicrun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
./bin/gr -i ./tests/aliceharmonicrun/SH_2pi_REF_ALICE.json -n 10000000
./bin/gr -i ./tests/aliceharmonicrun/SH_2pi_ALICE.json     -n 10000000

fi

# ***********************************************************************
# Fiducial cuts <ETAMIN,ETAMAX,PTMIN,PTMAX>
FIDCUTS=-0.9,0.9,0.1,100.0
# ***********************************************************************

# System kinematic variables binning <bins,min,max>
MBINS=30,0.25,1.75
PBINS=1,0.0,1.75
YBINS=1,-0.9,0.9

# PARAMETERS
LMAX=4
REMOVEODD=true
REMOVENEGATIVE=true
SVDREG=1e-4
L1REG=0 #1e-5
EML=true

# Lorentz frames
for FRAME in HE CS PG SR
do

# Analyze
#./bin/fitharmonic -r SH_2pi_REF_ALICE -i SH_2pi_ALICE -s true -t MC \
./bin/fitharmonic -r SH_2pi_REF_ALICE -i ALICE7_2pi.csv -t DATA -s false \
-c $FIDCUTS \
-f $FRAME -l $LMAX -o $REMOVEODD -n $REMOVENEGATIVE -a $SVDREG -b $L1REG -e $EML \
-M $MBINS -P $PBINS -Y $YBINS \
-s $FS \
-X 100000

done

# Implement 2D harmonic plots (M,Pt)
# ...
