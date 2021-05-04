#!/bin/bash
#
# Rivet analysis example of exclusive lepton pair production
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

EVENTS=50000

./bin/gr -i ./fitcard/STAR_1792394_pipi.json -n $EVENTS -w true -l true  \
-o "STAR_1792394_pipi_screened" -f "hepmc3"

fi


# Set RIVET environment here
source $HOME/local/rivet/rivetenv.sh

# Analyze and plot
SET=STAR_2020_I1792394

rivet ./output/STAR_1792394_pipi_screened.hepmc3 --analysis=$SET -o STAR_pipi_screened.yoda

rivet-mkhtml \
STAR_pipi_screened.yoda:"GRANIITTI (screened)" \
--mc-errs -o ./rivetplots/STAR_1792394_pipi/

