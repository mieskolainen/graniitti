#!/bin/bash
#
# Simple minimum bias simulation and RIVET analysis
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Not set
if [ -z "$EVENTS" ]
then
	EVENTS=10000
fi

./bin/minbias 900   $EVENTS
./bin/minbias 7000  $EVENTS
./bin/minbias 13000 $EVENTS

fi

# Set source here
source $HOME/scratch/rivet/local/rivetenv.sh


# # # # # -----------------------------------------------------------------------------------
# # # # # 0.96 TeV

# # Analyze and plot
# #SET=MC_PRINTEVENT
SET=ATLAS_2010_S8918562
#,ALICE_2010_S8625980,ALICE_2010_S8624100,ALICE_2010_S8706239,ALICE_2011_S8945144

# #--cross-section=7.36e10
rivet --analysis=$SET ./output/minbias_900.hepmc2
rivet-mkhtml --mc-errs Rivet.yoda:"GRANIITTI" -o ./rivetplots/minbias_900/


# # # # # -----------------------------------------------------------------------------------
# # # # # 7 TeV

# # # # # # Analyze and plot
SET=ATLAS_2010_S8918562,ALICE_2010_S8625980,ATLAS_2012_I1084540,TOTEM_2012_I1115294,CMS_2015_I1356998,ALICE_2015_I1357424
# # # # # #,CMS_2010_S8656010,ATLAS_2010_S8894728

#rivet --analysis=$SET ./output/minbias_7000.hepmc2
#rivet-mkhtml --mc-errs Rivet.yoda:"GRANIITTI" -o ./rivetplots/minbias_7000/


# # # # # -----------------------------------------------------------------------------------
# # # # # 13 TeV

# # # # # Analyze and plot
SET=ATLAS_2016_I1467230
#,CMS_2015_I1384119

rivet --analysis=$SET ./output/minbias_13000.hepmc2
rivet-mkhtml --mc-errs Rivet.yoda:"GRANIITTI" -o ./rivetplots/minbias_13000/
