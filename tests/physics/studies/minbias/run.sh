#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simple minimum bias simulation and RIVET analysis
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/minbias/run.sh

# Load an optional Rivet setup when it is not already active
study_require_opt_in RUN_RIVET "minimum bias Rivet study"
if [[ -n "${RIVET_ENV:-}" ]]
then
	source "${RIVET_ENV}"
fi
study_require_command rivet "minimum bias Rivet study"
study_require_command rivet-mkhtml "minimum bias Rivet study"

if study_generates
then

./bin/minbias 900   $NEVENTS
./bin/minbias 7000  $NEVENTS
./bin/minbias 13000 $NEVENTS

fi


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
