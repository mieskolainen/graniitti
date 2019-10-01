#!/bin/sh
#
# Loop over fiducial measurements and construct latex array
#
# Run with: source ./tests/exloop.sh



read -p "run: Screening on? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then
SCREENING=1
else
SCREENING=0
fi

echo ''
echo 'Hint: run cat /tmp/GRANIITTI.out or cat /tmp/GRANIITTI_S2.out for the progress ...'
echo ''
echo ''

echo '\hline \\'
echo 'Measurement' '&' 'xs +- stat +- syst' '&' 'screened xs' '&' 'unscreened xs' '&' '<S^2>' '\\'
echo '\hline \\'
echo ''

# Loop over steering cards
FILES=./tests/LHC_TEVATRON_RHIC/*.json
for F in $FILES
do

	# Extract filename + strip of extension
	CARD=$(basename $F)
	CARD="${CARD%.*}"

	MEASUREMENT=(0.0 0.0 0.0)
	
	# ====================================================================

	# http://cds.cern.ch/record/2679648/files/FSQ-16-006-pas.pdf
	if [ "$CARD" == "CMS19_2pi" ]
	then
	MEASUREMENT=(19.0E-6 0.6E-6 3.2E-6)
	fi
	#  <ALICE 7 TeV PWA |y(pi+pi-)| < 0.9
	if [ "$CARD" == "ALICE19_2pi_PWA" ]
	then
	MEASUREMENT=(31E-6 0.5E-6 2E-6)
	fi
	# https://discoverycenter.nbi.ku.dk/teaching/thesis_page/MasterEmilBolsFinal.pdf
	if [ "$CARD" == "ATLAS17_2pi" ]
	then
	MEASUREMENT=(18.75E-6 0.048E-6 0.770E-6)
	fi
	# https://arxiv.org/pdf/hep-ex/0611040.pdf
	if [ "$CARD" == "CDF07_ee" ]
	then
	MEASUREMENT=(1.6E-12 0.5E-12 0.3E-12)
	fi
	# https://arxiv.org/pdf/1112.0858.pdf
	if [ "$CARD" == "CDF11_ee" ]
	then
	MEASUREMENT=(2.88E-12 0.57E-12 0.63E-12)
	fi
	# https://arxiv.org/abs/1111.5536
	if [ "$CARD" == "CMS11_mumu" ]
	then
	MEASUREMENT=(3.38E-12 0.58E-12 0.21E-12)
	fi
	# https://arxiv.org/abs/1506.07098
	if [ "$CARD" == "ATLAS15_ee" ]
	then
	MEASUREMENT=(0.428E-12 0.035E-12 0.018E-12)
	fi
	# https://arxiv.org/abs/1506.07098
	if [ "$CARD" == "ATLAS15_mumu" ]
	then
	MEASUREMENT=(0.628E-12 0.032E-12 0.021E-12)
	fi
	# https://arxiv.org/abs/1708.04053
	if [ "$CARD" == "ATLAS17_mumu" ]
	then
	MEASUREMENT=(3.12E-12 0.07E-12 0.14E-12)
	fi

	# ====================================================================


	# Generate
	./bin/gr -i ./tests/LHC_TEVATRON_RHIC/$CARD.json -l false -n 0 > /tmp/GRANIITTI.out

	if [ $SCREENING == 1 ]
	then
	./bin/gr -i ./tests/LHC_TEVATRON_RHIC/$CARD.json -l true -n 0 > /tmp/GRANIITTI_S2.out		
	fi

	# Extract integrated cross section
	XS=$(grep -E 'barn' /tmp/GRANIITTI.out)
	XS=${XS#*[}
	XS=${XS%]*}
	XS=${XS%+-*} # Remove stat errors

	# Extract screened integrated cross section
	XS_S2=0.0
	if [ "$SCREENING" == 1 ]
	then
	XS_S2=$(grep -E 'barn' /tmp/GRANIITTI_S2.out)
	XS_S2=${XS_S2#*[}
	XS_S2=${XS_S2%]*}
	XS_S2=${XS_S2%+-*} # Remove stat errors
	fi

	# Add +- characters
	printf -v ids_d ' +- %s' "${MEASUREMENT[@]}" # yields "value +- stat +- syst"
	ids_d=${ids_d:3}                             # remove the leading ' +- '


	# Calculate division
	S2=$(python -c "print('%0.2f' % ($XS_S2/$XS))")


	# Print out
	echo $CARD ' & ' $ids_d ' & ' $XS_S2 ' & ' $XS ' & ' $S2 '\\'
	#printf '%s & %s & %s & %s \n' $CARD $ids_d $XS_S2 $XS

done

echo ''
echo 'exloop.sh [DONE]'


# Analyze
#  // <https://indico.cern.ch/event/713101/contributions/3102315/attachments/1705771/2748440/Diffraction2018_RafalSikora.pdf>
#  input.push_back(MEASUREMENT(PATH + "STAR18_2pi.json" 0 0 0));

