#!/bin/bash
#
# Rivet analysis example of exlusive lepton pair production
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate events
./bin/gr ./tests/LHC_TEVATRON_RHIC/CMS11_mumu.json

fi

# Set source here (otherwise collapse between HepMC2 and HepMC3 libraries!)

#source ./lenovo/rivet/local/rivetenv.sh
source ./i5/local/rivetenv.sh

# -----------------------------------------------------------------------------------
# 7 TeV

# Analyze and plot
SET=CMS_2011_I954992

rivet --analysis=$SET ./output/CMS11_mumu.hepmc2
rivet-mkhtml --mc-errs Rivet.yoda:"GRANIITTI" -o ./rivetplots/CMS11_mumu/

