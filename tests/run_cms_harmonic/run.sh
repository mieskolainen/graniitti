#!/bin/sh
#
# Simulation pi+pi- and spherical harmonic expansion with CMS cuts
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Not set
if [ -z "$EVENTS" ]
then
	EVENTS=1000000
fi

# Generate
./bin/gr -i ./tests/run_cms_harmonic/SH_2pi_J0_CMS.json -n $EVENTS
./bin/gr -i ./tests/run_cms_harmonic/SH_2pi_CMS.json    -n $EVENTS

fi

# ***********************************************************************
# Fiducial cuts <ETAMIN,ETAMAX,PTMIN,PTMAX> (no white space!)
FIDCUTS=-2.5,2.5,0.1,100.0
# ***********************************************************************

# System kinematic variables binning <bins,min,max> (no white space!)
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
for FRAME in CM HX CS PG GJ
do

# Expand the data
./bin/fitharmonic \
-r "SH_2pi_J0_CMS" \
-i "SH_2pi_J0_CMS, SH_2pi_CMS" \
-l "GRANIITTI J=0, GRANIITTI #pi^{+}#pi^{-}" \
-d "MC, MC" \
-z "true, true" \
-t "#Omega{Detector}: |#eta| < 2.5 #wedge p_{T} > 0.1 GeV,#Omega{Fiducial}: |#eta| < 2.5 #wedge p_{T} > 0.1 GeV,#Omega{Flat}: |Y_{x}| < 2.5" \
-c $FIDCUTS \
-f $FRAME -g $LMAX -o $REMOVEODD -v $REMOVENEGATIVE -a $SVDREG -b $L1REG -e $EML \
-M $MBINS -P $PBINS -Y $YBINS \
-S "-1.0, -1.0" \
-w "" \
-X 1000000000

# -S -1 for normalized rates

done

# Implement 2D harmonic plots (M,Pt)
# ...