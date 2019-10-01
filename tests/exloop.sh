#!/bin/sh
#
# Loop over measurements
#
# Run with: source ./tests/exloop.sh

read -p "run: Screening on? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then
SCREENING="true"
else
SCREENING="false"
fi

# Loop over steering cards
FILES=./tests/LHC_TEVATRON_RHIC/*.json
for F in $FILES
do

	# Extract filename + strip of extension
	CARD=$(basename $F)
	CARD="${CARD%.*}"

	MEASUREMENT=(0.0, 0.0, 0.0)
	
	# http://cds.cern.ch/record/2679648/files/FSQ-16-006-pas.pdf
	if [ "$CARD" == "CMS19_2pi" ]
	then
	MEASUREMENT=(19.0e-6, 0.6e-6, 3.2e-6)
	fi
	#  <ALICE 7 TeV PWA, |y(pi+pi-)| < 0.9
	if [ "$CARD" == "ALICE19_2pi_PWA" ]
	then
	MEASUREMENT=(31e-6, 0.5e-6, 2e-6)
	fi
	# https://discoverycenter.nbi.ku.dk/teaching/thesis_page/MasterEmilBolsFinal.pdf
	if [ "$CARD" == "ATLAS17_2pi" ]
	then
	MEASUREMENT=(18.75e-6, 0.048e-6, 0.770e-6)
	fi
	# https://arxiv.org/pdf/hep-ex/0611040.pdf
	if [ "$CARD" == "CDF07_ee" ]
	then
	MEASUREMENT=(1.6e-12, 0.5e-12, 0.3e-12)
	fi
	# https://arxiv.org/pdf/1112.0858.pdf
	if [ "$CARD" == "CDF11_ee" ]
	then
	MEASUREMENT=(2.88e-12, 0.57e-12, 0.63e-12)
	fi
	# https://arxiv.org/abs/1111.5536
	if [ "$CARD" == "CMS11_mumu" ]
	then
	MEASUREMENT=(3.38e-12, 0.58e-12, 0.21e-12)
	fi
	# https://arxiv.org/abs/1506.07098
	if [ "$CARD" == "ATLAS15_ee" ]
	then
	MEASUREMENT=(0.428e-12, 0.035e-12, 0.018e-12)
	fi
	# https://arxiv.org/abs/1506.07098
	if [ "$CARD" == "ATLAS15_mumu" ]
	then
	MEASUREMENT=(0.628e-12, 0.032e-12, 0.021e-12)
	fi
	# https://arxiv.org/abs/1708.04053
	if [ "$CARD" == "ATLAS17_mumu" ]
	then
	MEASUREMENT=(3.12e-12, 0.07e-12, 0.14e-12)
	fi

	# Generate
	echo '['$CARD']'
	./bin/gr -i ./tests/LHC_TEVATRON_RHIC/$CARD.json -l $SCREENING -n 0 > /tmp/GRANIITTI.out

	SIMULATION=$(grep -E 'barn' /tmp/GRANIITTI.out)
	echo $SIMULATION ',' ${MEASUREMENT[@]}
	echo ''

done

# Analyze
#  // <https://indico.cern.ch/event/713101/contributions/3102315/attachments/1705771/2748440/Diffraction2018_RafalSikora.pdf>
#  input.push_back(MEASUREMENT(PATH + "STAR18_2pi.json", 0, 0, 0));

#  // <ATLAS W+W->
#  input.push_back(MEASUREMENT(PATH + "ATLAS_WW.json", 0, 0, 0));
