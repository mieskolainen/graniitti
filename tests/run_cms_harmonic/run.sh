#!/bin/sh
#
# Simulation and spherical harmonic expansion with CMS cuts
#
# Run with: source ./tests/run_xxx/run.sh

read -p "cmsharmonicrun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
./bin/gr -i ./tests/run_cms_harmonic/SH_2pi_J0_CMS.json -n 1000000
./bin/gr -i ./tests/run_cms_harmonic/SH_2pi_CMS.json    -n 1000000

fi

# ***********************************************************************
# Fiducial cuts <ETAMIN,ETAMAX,PTMIN,PTMAX>
FIDCUTS=-2.5,2.5,0.1,100.0
# ***********************************************************************

# System kinematic variables binning <bins,min,max>
MBINS=30,0.28,2.0
PBINS=1,0.0,10.0
YBINS=1,-2.5,2.5

# PARAMETERS
LMAX=4
REMOVEODD=true
REMOVENEGATIVE=true
SVDREG=1e-5
L1REG=0
EML=false

# Lorentz frames
for FRAME in HE CS PG SR
do

# Expand the data
./bin/fitharmonic -r SH_2pi_J0_CMS -i SH_2pi_J0_CMS,SH_2pi_CMS \
-l 'GRANIITTI J=0,GRANIITTI #pi^{+}#pi^{-}' \
-d MC,MC \
-z true,true \
-t '#Omega{Detector}: |#eta| < 2.5 #wedge p_{T} > 0.1 GeV,#Omega{Fiducial}: |#eta| < 2.5 #wedge p_{T} > 0.1 GeV,#Omega{Flat}: |Y_{x}| < 2.5' \
-c $FIDCUTS \
-f $FRAME -g $LMAX -o $REMOVEODD -v $REMOVENEGATIVE -a $SVDREG -b $L1REG -e $EML \
-M $MBINS -P $PBINS -Y $YBINS \
-S -1.0,-1.0 \
-X 100000

# -S -1 for normalized rates

done

# Implement 2D harmonic plots (M,Pt)
# ...