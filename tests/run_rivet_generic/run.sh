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
./bin/gr -i ./tests/LHC_TEVATRON_RHIC/CMS11_mumu.json -n 1000 -w false -l false -o "CMS11_mumu" -f "hepmc3"

fi

# Set source here
source $HOME/local/rivet/rivetenv.sh


# Analyze and plot
SET=MC_GENERIC

rivet ./output/CMS11_mumu.hepmc3 --analysis=$SET
rivet-mkhtml --mc-errs Rivet.yoda:"GRANIITTI" -o ./rivetplots/MC_GENERIC/


#usage: rivet [-h] [--version] [-a ANA] [--list-analyses] [--list-keywords] [--list-used-analyses] [--show-analysis SHOW_ANALYSES] [--show-bibtex] [--analysis-path PATH] [--analysis-path-append PATH]
#             [--pwd] [-o HISTOFILE] [-p PRELOADFILE] [--no-histo-file] [-x XS] [-n NUM] [--nskip NUM] [--skip-weights] [--match-weights MATCH_WEIGHTS] [--unmatch-weights UNMATCH_WEIGHTS]
#             [--nominal-weight NOMINAL_WEIGHT] [--weight-cap WEIGHT_CAP] [--nlo-smearing NLO_SMEARING] [--runname NAME] [--ignore-beams] [-d NUM] [--event-timeout NSECS] [--run-timeout NSECS]
#             [-l NATIVE_LOG_STRS] [-v] [-q]
#             [ARGS [ARGS ...]]

