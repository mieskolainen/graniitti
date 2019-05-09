#!/bin/bash
#
# Simple minimum bias simulation and RIVET analysis
#
# Run with: source ./tests/run_xxx/run.sh

./bin/minbias 900   10000
./bin/minbias 7000  10000
./bin/minbias 13000 10000


# SET RIVET source here (otherwise collapse between HepMC2 and HepMC3 libraries!)
source ./lenovo/rivet/local/rivetenv.sh
#source ./i5/local/rivetenv.sh
#source ./core4/rivet/local/rivetenv.sh


# # # # # -----------------------------------------------------------------------------------
# # # # # 0.96 TeV

# # Analyze and plot
# #SET=MC_PRINTEVENT
SET=ATLAS_2010_S8918562
#,ALICE_2010_S8625980,ALICE_2010_S8624100,ALICE_2010_S8706239,ALICE_2011_S8945144

# #--cross-section=7.36e10
rivet --analysis=$SET ./output/minbias_900.hepmc2
rivet-mkhtml --mc-errs Rivet.yoda:"GRANIITTI" -o ./rivet-plots-900/


# # # # # -----------------------------------------------------------------------------------
# # # # # 7 TeV

# # # # # # Analyze and plot
SET=ATLAS_2010_S8918562,ALICE_2010_S8625980,ATLAS_2012_I1084540,TOTEM_2012_I1115294,CMS_2015_I1356998,ALICE_2015_I1357424
# # # # # #,CMS_2010_S8656010,ATLAS_2010_S8894728

#rivet --analysis=$SET ./output/minbias_7000.hepmc2
#rivet-mkhtml --mc-errs Rivet.yoda:"GRANIITTI" -o ./rivet-plots-7000/


# # # # # -----------------------------------------------------------------------------------
# # # # # 13 TeV

# # # # # Analyze and plot
SET=ATLAS_2016_I1467230
#,CMS_2015_I1384119

rivet --analysis=$SET ./output/minbias_13000.hepmc2
rivet-mkhtml --mc-errs Rivet.yoda:"GRANIITTI" -o ./rivet-plots-13000/
